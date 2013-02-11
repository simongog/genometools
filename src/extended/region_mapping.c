/*
  Copyright (c) 2007-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2012-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2012 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "core/assert_api.h"
#include "core/bioseq.h"
#include "core/bioseq_col.h"
#include "core/encseq_col.h"
#include "core/ma.h"
#include "core/md5_seqid.h"
#include "core/seq_col.h"
#include "core/str_array.h"
#include "core/undef_api.h"
#include "extended/mapping.h"
#include "extended/region_mapping_api.h"
#include "extended/seqid2seqnum_mapping.h"

struct GtRegionMapping {
  GtStrArray *sequence_filenames;
  GtStr *sequence_file,  /* the (current) sequence file */
        *sequence_name;  /* the (current) sequence name */
  bool matchdesc,
       usedesc,
       userawseq;
  GtMapping *mapping;
  GtBioseq *bioseq; /* the current bioseq */
  GtEncseq *encseq;
  GtSeqCol *seq_col;
  GtSeqid2SeqnumMapping *seqid2seqnum_mapping;
  const char *rawseq;
  unsigned long rawlength,
                rawoffset;
  unsigned int reference_count;
};

GtRegionMapping* gt_region_mapping_new_mapping(GtStr *mapping_filename,
                                               GtError *err)
{
  GtRegionMapping *rm;
  gt_error_check(err);
  gt_assert(mapping_filename);
  rm = gt_calloc(1, sizeof (GtRegionMapping));
  rm->mapping = gt_mapping_new(mapping_filename, "mapping",
                               GT_MAPPINGTYPE_STRING, err);
  if (!rm->mapping) {
    gt_region_mapping_delete(rm);
    return NULL;
  }
  return rm;
}

GtRegionMapping* gt_region_mapping_new_seqfiles(GtStrArray *sequence_filenames,
                                                bool matchdesc, bool usedesc)
{
  GtRegionMapping *rm;
  gt_assert(sequence_filenames);
  gt_assert(!(matchdesc && usedesc));
  rm = gt_calloc(1, sizeof (GtRegionMapping));
  rm->sequence_filenames = gt_str_array_ref(sequence_filenames);
  rm->matchdesc = matchdesc;
  rm->usedesc = usedesc;
  return rm;
}

GtRegionMapping* gt_region_mapping_new_encseq(GtEncseq *encseq, bool matchdesc,
                                              bool usedesc)
{
  GtRegionMapping *rm;
  gt_assert(encseq);
  gt_assert(!(matchdesc && usedesc));
  rm = gt_calloc(1, sizeof (GtRegionMapping));
  rm->encseq = gt_encseq_ref(encseq);
  rm->matchdesc = matchdesc;
  rm->usedesc = usedesc;
  return rm;
}

GtRegionMapping* gt_region_mapping_new_rawseq(const char *rawseq,
                                              unsigned long length,
                                              unsigned long offset)
{
  GtRegionMapping *rm;
  gt_assert(rawseq);
  rm = gt_calloc(1, sizeof (GtRegionMapping));
  rm->userawseq = true;
  rm->rawseq = rawseq;
  rm->rawlength = length;
  rm->rawoffset = offset;
  return rm;
}

GtRegionMapping* gt_region_mapping_ref(GtRegionMapping *rm)
{
  gt_assert(rm);
  rm->reference_count++;
  return rm;
}

static GtStr* region_mapping_map(GtRegionMapping *rm,
                                 const char *sequence_region, GtError *err)
{
  gt_error_check(err);
  gt_assert(rm && sequence_region);
  if (!rm->mapping)
    return gt_str_ref(gt_str_array_get_str(rm->sequence_filenames, 0));
  else
    return gt_mapping_map_string(rm->mapping, sequence_region, err);
}

static int update_seq_col_if_necessary(GtRegionMapping *rm, GtStr *seqid,
                                       GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(rm && seqid);
  /* for mappings, we need to load the changed sequence, if needed... */
  if (rm->mapping) {
    if (!rm->sequence_file || (gt_str_cmp(rm->sequence_name, seqid))) {
      gt_str_delete(rm->sequence_file);
      /* ignore MD5 hashes when using region mappings */
      if (gt_md5_seqid_has_prefix(gt_str_get(seqid))) {
        rm->sequence_file = region_mapping_map(rm,
                                               gt_str_get(seqid)
                                                 +GT_MD5_SEQID_TOTAL_LEN,
                                               err);
      } else
        rm->sequence_file = region_mapping_map(rm, gt_str_get(seqid), err);
      if (!rm->sequence_file)
        had_err = -1;
      else {
        /* load new seqcol */
        if (!rm->sequence_filenames)
          rm->sequence_filenames = gt_str_array_new();
        else
          gt_str_array_reset(rm->sequence_filenames);
        gt_str_array_add(rm->sequence_filenames, rm->sequence_file);
        if (!rm->sequence_name)
          rm->sequence_name = gt_str_new();
        else
          gt_str_reset(rm->sequence_name);
        gt_str_append_str(rm->sequence_name, seqid);
        gt_seq_col_delete(rm->seq_col);
        rm->seq_col = gt_bioseq_col_new(rm->sequence_filenames, err);
        if (!rm->seq_col)
          had_err = -1;
      }
    }
  } else {
    /* ...otherwise, just make sure the seqcol is loaded */
    if (!rm->seq_col) {
      if (rm->encseq) {
        if (!(rm->seq_col = gt_encseq_col_new(rm->encseq, err)))
          had_err = -1;
      } else {
        gt_assert(rm->sequence_filenames);
        if (!(rm->seq_col = gt_bioseq_col_new(rm->sequence_filenames, err)))
          had_err = -1;
      }
    }
    if (!had_err && rm->usedesc) {
      if (rm->seqid2seqnum_mapping)
        gt_seqid2seqnum_mapping_delete(rm->seqid2seqnum_mapping);
      rm->seqid2seqnum_mapping =
                           gt_seqid2seqnum_mapping_new_seqcol(rm->seq_col, err);
      if (!rm->seqid2seqnum_mapping) {
        had_err = -1;
      }
    }
  }
  return had_err;
}

int gt_region_mapping_get_sequence(GtRegionMapping *rm, char **seq,
                                   GtStr *seqid, unsigned long start,
                                   unsigned long end, GtError *err)
{
  int had_err = 0;
  unsigned long offset = 1;
  GtRange range = {GT_UNDEF_ULONG, GT_UNDEF_ULONG};
  gt_error_check(err);
  gt_assert(rm && seq && seqid && gt_str_length(seqid) > 0);

  /* handle rawseq access first  */
  if (rm->userawseq) {
    gt_assert(!rm->seqid2seqnum_mapping);
    *seq = gt_calloc(end - start + 1, sizeof (char));
    strncpy(*seq, rm->rawseq + start - 1, (end - start + 1) * sizeof (char));
    return 0;
  }

  /* make sure that correct sequence is loaded */
  had_err = update_seq_col_if_necessary(rm, seqid, err);

  /* MD5 sequence id */
  if (!had_err) {
    if (gt_md5_seqid_has_prefix(gt_str_get(seqid))) {
      had_err = gt_seq_col_md5_to_seq(rm->seq_col, seq, start - offset,
                                        end - offset, seqid, err);
      return had_err;
    }
  }

  /* ``regular'' sequence ID */
  if (!had_err) {
    gt_assert(!rm->usedesc || rm->seqid2seqnum_mapping);
    gt_assert(rm->mapping || rm->seq_col);
    if (rm->usedesc) {
      unsigned long seqnum, filenum;
      gt_assert(rm->seqid2seqnum_mapping);
      range.start = start;
      range.end = end;
      had_err = gt_seqid2seqnum_mapping_map(rm->seqid2seqnum_mapping,
                                            gt_str_get(seqid), &range, &seqnum,
                                            &filenum, &offset, err);

      if (!had_err) {
        if (range.end != GT_UNDEF_ULONG && range.start != GT_UNDEF_ULONG &&
              range.end >= gt_seq_col_get_sequence_length(rm->seq_col, filenum,
                                                          seqnum)
              + offset) {
          gt_error_set(err, "trying to extract range %lu-%lu on sequence "
                       "``%s''which is not covered by that sequence (with "
                       "boundaries %lu-%lu). Has the sequence-region "
                       "to sequence mapping been defined correctly?",
                       start, end, gt_str_get(seqid),
                       range.start, range.end);
          had_err = -1;
        }
      }
      if (!had_err) {
        *seq = gt_seq_col_get_sequence(rm->seq_col, filenum, seqnum,
                                       start - offset, end - offset);
      }
    } else if (rm->matchdesc) {
      gt_assert(!rm->seqid2seqnum_mapping);
      gt_assert(rm->seq_col);
      if (!had_err) {
        had_err = gt_seq_col_grep_desc(rm->seq_col, seq, start - 1, end - 1,
                                       seqid, err);
      }
    } else if (rm->userawseq) {
      gt_assert(!rm->seqid2seqnum_mapping);
      *seq = gt_calloc(end - start + 1, sizeof (char));
      strncpy(*seq, rm->rawseq + start - 1, (end - start + 1) * sizeof (char));
    } else {
      gt_assert(rm->seq_col);

      if (!had_err) {
        unsigned long seqlength = gt_seq_col_get_sequence_length(rm->seq_col,
                                                                 0, 0);
        if (start > seqlength || end > seqlength) {
          had_err = -1;
          gt_error_set(err, "trying to extract range %lu-%lu on sequence "
                       "``%s''which is not covered by that sequence (only "
                       "%lu characters in size). Has the sequence-region "
                       "to sequence mapping been defined correctly?",
                       start, end, gt_str_get(seqid), seqlength);
        }
        if (!had_err) {
          *seq = gt_seq_col_get_sequence(rm->seq_col, 0, 0, start - offset,
                                         end - offset);
        }
      }
    }
  }
  return had_err;
}

int gt_region_mapping_get_sequence_length(GtRegionMapping *rm,
                                          unsigned long *length, GtStr *seqid,
                                          GtError *err)
{
  int had_err;
  unsigned long filenum, seqnum, GT_UNUSED offset;
  gt_error_check(err);
  GT_UNUSED GtRange range;
  gt_assert(rm && seqid);
  if (rm->userawseq) {
    return rm->rawlength;
  }
  had_err = update_seq_col_if_necessary(rm, seqid, err);
  if (!had_err) {
    if (gt_md5_seqid_has_prefix(gt_str_get(seqid))) {
      had_err = gt_seq_col_md5_to_sequence_length(rm->seq_col, length, seqid,
                                                  err);
    }
    return had_err;
  }
  if (!had_err) {
    if (rm->usedesc) {
      gt_assert(rm->seqid2seqnum_mapping);
      had_err = gt_seqid2seqnum_mapping_map(rm->seqid2seqnum_mapping,
                                            gt_str_get(seqid), &range, &seqnum,
                                            &filenum,
                                            &offset, err);
      if (!had_err)
        *length = gt_seq_col_get_sequence_length(rm->seq_col, filenum, seqnum);
    }
    else if (rm->matchdesc) {
      had_err = gt_seq_col_grep_desc_sequence_length(rm->seq_col, length,
                                                     seqid, err);
    }
    else {
      if (!had_err)
        *length = gt_seq_col_get_sequence_length(rm->seq_col, 0, 0);
    }
  }
  return had_err;
}

int gt_region_mapping_get_description(GtRegionMapping *rm, GtStr *desc,
                                      GtStr *md5_seqid, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(rm && desc && md5_seqid);
  /* this method is only implemented for MD5 seqids */
  gt_assert(gt_md5_seqid_has_prefix(gt_str_get(md5_seqid)));
  had_err = update_seq_col_if_necessary(rm, md5_seqid, err);
  if (!had_err)
    had_err = gt_seq_col_md5_to_description(rm->seq_col, desc, md5_seqid, err);
  return had_err;
}

const char* gt_region_mapping_get_md5_fingerprint(GtRegionMapping *rm,
                                                  GtStr *seqid,
                                                  const GtRange *range,
                                                  unsigned long *offset,
                                                  GtError *err)
{
  const char *md5 = NULL;
  int had_err;
  unsigned long filenum, seqnum;
  gt_error_check(err);
  gt_assert(rm && seqid);
  gt_assert(!rm->userawseq); /* not implemented */
  gt_assert(!gt_md5_seqid_has_prefix(gt_str_get(seqid))); /* not implemented */
  had_err = update_seq_col_if_necessary(rm, seqid, err);
  if (!had_err) {
    if (rm->usedesc) {
      gt_assert(rm->seqid2seqnum_mapping);
      had_err = gt_seqid2seqnum_mapping_map(rm->seqid2seqnum_mapping,
                                            gt_str_get(seqid), range, &seqnum,
                                            &filenum,
                                            offset, err);
      if (!had_err)
        md5 = gt_seq_col_get_md5_fingerprint(rm->seq_col, filenum, seqnum);
    }
    else if (rm->matchdesc) {
      if (!rm->seq_col) {
        if (rm->encseq) {
          if (!(rm->seq_col = gt_encseq_col_new(rm->encseq, err)))
            had_err = -1;
        } else {
          if (!(rm->seq_col = gt_bioseq_col_new(rm->sequence_filenames, err)))
            had_err = -1;
        }
      }
      if (!had_err)
        had_err = gt_seq_col_grep_desc_md5(rm->seq_col, &md5, seqid, err);
      *offset = 1;
    }
    else {
      if (!had_err)
        md5 = gt_seq_col_get_md5_fingerprint(rm->seq_col, 0, 0);
      *offset = 1;
    }
  }
  return md5;
}

void gt_region_mapping_delete(GtRegionMapping *rm)
{
  if (!rm) return;
  if (rm->reference_count) {
    rm->reference_count--;
    return;
  }
  gt_str_array_delete(rm->sequence_filenames);
  gt_str_delete(rm->sequence_file);
  gt_str_delete(rm->sequence_name);
  gt_mapping_delete(rm->mapping);
  gt_bioseq_delete(rm->bioseq);
  gt_encseq_delete(rm->encseq);
  gt_seq_col_delete(rm->seq_col);
  gt_seqid2seqnum_mapping_delete(rm->seqid2seqnum_mapping);
  gt_free(rm);
}
