/*
  Copyright (c) 2012      Manuela Beckert <9beckert@informatik.uni-hamburg.de>
  Copyright (c) 2012      Dorle Osterode <9osterod@informatik.uni-hamburg.de>
  Copyright (c) 2012-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012-2013 Center for Bioinformatics, University of Hamburg

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

#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "core/arraydef.h"
#include "core/encseq.h"
#include "core/log_api.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/str_api.h"
#include "core/undef_api.h"
#include "extended/genome_node.h"
#include "extended/node_stream_api.h"
#include "extended/region_node_api.h"
#include "extended/tir_stream.h"
#include "match/esa-maxpairs.h"
#include "match/esa-mmsearch.h"
#include "match/esa-seqread.h"
#include "match/greedyedist.h"
#include "match/querymatch.h"
#include "match/xdrop.h"
#include "extended/tir_stream.h"

typedef struct
{
  unsigned long pos1;         /* position of first seed */
  unsigned long pos2;         /* position of second seed (other contig) */
  unsigned long offset;       /* distance between them related to the actual
                                 sequence (not mirrored) */
  unsigned long len;          /* length of the seed  */
  unsigned long contignumber; /* number of contig for this seed */
} Seed;
GT_DECLAREARRAYSTRUCT(Seed);

typedef struct
{
  GtArraySeed seed;
  unsigned long max_tir_length;
  unsigned long min_tir_length;
  unsigned long max_tir_distance;
  unsigned long min_tir_distance;
  unsigned long num_of_contigs;
  unsigned long midpos;
  unsigned long totallength;
} SeedInfo;

typedef struct
{
  unsigned long contignumber,
                left_tir_start,  /* first position of TIR on forward strand */
                left_tir_end,    /* last position of TIR on forward strand */
                right_tir_start, /* first position of TIR on reverse strand */
                right_tir_end;   /* last position of TIR on reverse strand */
  double        similarity;      /* similarity of the two TIRs */
  bool          skip;            /* needed to remove overlaps if wanted */
  unsigned long tsd_length;      /* length of tsd at start of left tir and end
                                    of right tir */
} TIRPair;

GT_DECLAREARRAYSTRUCT(TIRPair);

typedef struct
{
  unsigned long left_start_pos,   /* represents the start position for the TSD
                                     search */
                right_start_pos;
  GtArraySeed TSDs;               /* array to store the TSDs */
} TSDinfo;

typedef enum {
  GT_TIR_STREAM_STATE_START,
  GT_TIR_STREAM_STATE_REGIONS,
  GT_TIR_STREAM_STATE_COMMENTS,
  GT_TIR_STREAM_STATE_FEATURES
} GtTIRStreamState;

struct GtTIRStream
{
  const GtNodeStream          parent_instance;
  const GtEncseq              *encseq;
  Sequentialsuffixarrayreader *ssar;
  SeedInfo                    seedinfo;
  GtArrayTIRPair              tir_pairs;
  GtTIRStreamState            state;

  unsigned long               num_of_tirs,
                              cur_elem_index,
                              prev_seqnum;

  /* options */
  GtStr                       *str_indexname;
  unsigned long               min_seed_length,
                              min_TIR_length,
                              max_TIR_length,
                              min_TIR_distance,
                              max_TIR_distance;
  GtXdropArbitraryscores      arbit_scores;
  int                         xdrop_belowscore;
  double                      similarity_threshold;
  bool                        no_overlaps;
  bool                        best_overlaps;
  unsigned long               min_TSD_length,
                              max_TSD_length,
                              vicinity;
};

static int gt_tir_store_seeds(void *info, const GtEncseq *encseq,
                              unsigned long len, unsigned long pos1,
                              unsigned long pos2, GtError *err)
{
  unsigned long seqnum1,
                seqnum2;
  Seed *nextfreeseedpointer;
  unsigned long distance;
  SeedInfo *seeds = (SeedInfo *) info;
  gt_error_check(err);

  /* ensure order of seeds */
  if (pos1 > pos2) {
    unsigned long tmp = 0;
    tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }

  /* match mirrored vs. unmirrored */
  if (pos1 > seeds->midpos || pos2 < seeds->midpos)
    return 0;

  /* check distance constraints */
  distance = (GT_REVERSEPOS(seeds->totallength, pos2) - len + 1) - pos1;
  if (distance < seeds->min_tir_distance || distance > seeds->max_tir_distance)
    return 0;

  /* check whether matches are on the `same' contig */
  seqnum1 = gt_encseq_seqnum(encseq, pos1);
  seqnum2 = gt_encseq_seqnum(encseq, pos2);
  if (seqnum2 != seeds->num_of_contigs - seqnum1 - 1)
    return 0;

  /*char seq1[BUFSIZ], seq2[BUFSIZ];

  printf("ctg %lu | pos1 %lu | pos2 %lu | offset %lu | len %lu | midpos %lu\n",
         seqnum1, pos1, GT_REVERSEPOS(seeds->totallength, pos2) - len + 1,
         distance, len, seeds->midpos);
  gt_encseq_extract_decoded(encseq, seq1, pos1, pos1+len);
  seq1[len] = '\0';
  gt_encseq_extract_decoded(encseq, seq2,
                            GT_REVERSEPOS(seeds->totallength, pos2) - len + 1,
                            GT_REVERSEPOS(seeds->totallength, pos2) );
  seq2[len] = '\0';
  printf("%s  %s\n",seq1, seq2); */
  GT_GETNEXTFREEINARRAY(nextfreeseedpointer, &seeds->seed, Seed, 64);
  nextfreeseedpointer->pos1 = pos1;
  nextfreeseedpointer->pos2 = pos2;
  nextfreeseedpointer->offset = distance;
  nextfreeseedpointer->len = len;
  nextfreeseedpointer->contignumber = seqnum1;
  return 0;
}

/* this function is a call back function to store all TSDs found */
static int gt_tir_store_TSDs(void *info, GT_UNUSED const GtEncseq *encseq,
                             const GtQuerymatch *querymatch,
                             GT_UNUSED const GtUchar *query,
                             GT_UNUSED unsigned long query_totallength,
                             GT_UNUSED GtError *err)
{
  Seed *nextfree;
  TSDinfo *TSDs = (TSDinfo *) info;

  /* store the TSD at the next free index of info */
  GT_GETNEXTFREEINARRAY(nextfree, &TSDs->TSDs, Seed, 10);
  nextfree->pos1 = TSDs->left_start_pos + gt_querymatch_dbstart(querymatch);
  nextfree->offset = TSDs->right_start_pos
                       + gt_querymatch_querystart(querymatch)
                       - (nextfree->pos1);
  nextfree->len = gt_querymatch_querylen(querymatch);

  return 0;
}

static int gt_tir_compare_TIRs(TIRPair *pair1, TIRPair *pair2)
{
  if (pair1->contignumber < pair2->contignumber) {
    return -1;
  } else if (pair1->contignumber == pair2->contignumber) {
    if (pair1->left_tir_start < pair2->left_tir_start) {
      return -1;
    } else if (pair1->left_tir_start == pair2-> left_tir_start) {
      if (pair1->right_tir_start < pair2->right_tir_start) {
        return -1;
      } else if (pair1->right_tir_start == pair2->right_tir_start) {
        return 0;
      }
    }
  }
  return 1;
}

static void gt_tir_swap_TIRs(GtArrayTIRPair *tir_pairs, unsigned long pos1,
                             unsigned long pos2)
{
  TIRPair tmp;
  tmp = tir_pairs->spaceTIRPair[pos1];
  tir_pairs->spaceTIRPair[pos1] = tir_pairs->spaceTIRPair[pos2];
  tir_pairs->spaceTIRPair[pos2] = tmp;
}

/* sorts an array of tirs with inplace quicksort*/
static void gt_tir_sort_TIRs(GtArrayTIRPair *tir_pairs, unsigned long start,
                             unsigned long end)
{
  TIRPair pivot;
  unsigned long l;
  unsigned long r;

  if (end > start) {
    pivot = tir_pairs->spaceTIRPair[start];
    l = start + 1;
    r = end;
    while (l < r) {
      if (gt_tir_compare_TIRs(&tir_pairs->spaceTIRPair[l], &pivot) <= 0) {
        l++;
      } else {
        r--;
        gt_tir_swap_TIRs(tir_pairs, l, r);
      }
    }
    l--;
    gt_tir_swap_TIRs(tir_pairs, start, l);
    gt_tir_sort_TIRs(tir_pairs, start, l);
    gt_tir_sort_TIRs(tir_pairs, r, end);
  }
}

static unsigned long gt_tir_remove_overlaps(GtArrayTIRPair *src,
                                            GtArrayTIRPair *dest,
                                            unsigned long size_of_array,
                                            bool remove_all)
{
  TIRPair *new_pair,
          *pair1,
          *pair2;
  unsigned long start_pair1,
                end_pair1,
                start_pair2,
                end_pair2,
                num_of_tirs = size_of_array;
  int i, j;

  for (i = 0; i < size_of_array; i++) {
    pair1 = &src->spaceTIRPair[i];

    /* to prevent that a skipped one is checked twice */
    if (pair1->skip) continue;

    start_pair1 = pair1->left_tir_start;
    end_pair1 = pair1->right_tir_end;

    for (j = i + 1; j < size_of_array; j++) {
      pair2 = &src->spaceTIRPair[j];

      /* to prevent that a skipped one is checked twice */
      if (pair2->skip) continue;

      start_pair2 = pair2->left_tir_start;
      end_pair2 = pair2->right_tir_end;

      /* check if there is an overlap */
      if (!((end_pair1 < start_pair2) || (end_pair2 < start_pair1))) {
        /* no overlaps allowed */
        if (remove_all) {
          /* set skip to delete the tirs later */
          if (!pair1->skip) {
            num_of_tirs = num_of_tirs - 2;
          } else {
            num_of_tirs = num_of_tirs - 1;
          }
          pair1->skip = true;
          pair2->skip = true;
        } else {
          /* take the tir with best similarity */
          if (pair1->similarity >= pair2->similarity) {
            pair2->skip = true;
            num_of_tirs = num_of_tirs - 1;
          } else {
            pair1->skip = true;
            num_of_tirs = num_of_tirs - 1;
            break;
          }
        }
      }
    }
  }

  for (i = 0; i < size_of_array; i++) {
    pair1 = &src->spaceTIRPair[i];

    if (!pair1->skip) {
      GT_GETNEXTFREEINARRAY(new_pair, dest, TIRPair, 5);
      *new_pair = *pair1;
    }
  }

  return num_of_tirs;
}

static void gt_tir_find_best_TSD(TSDinfo *info, GtTIRStream *tir_stream,
                                 TIRPair *tir_pair)
{
  int i, j;
  Seed *tsd;
  unsigned long tsd_length,
                optimal_tsd_length,
                new_left_tir_start = tir_pair->left_tir_start,
                new_right_tir_end = tir_pair->right_tir_end,
                new_cost_left = 0,
                new_cost_right = 0,
                best_cost = ULONG_MAX,
                new_cost = 0;

  for (i = 0; i < info->TSDs.nextfreeSeed; i++) {

    tsd = &info->TSDs.spaceSeed[i];
    optimal_tsd_length = tsd->len;

    if (tsd->len <= tir_stream->min_TSD_length) continue;

    for (j = 0; j < tsd->len - tir_stream->min_TSD_length + 1; j++) {
      tsd_length = tsd->len - j;

      if (tsd_length < tir_stream->max_TSD_length) {
        if (tsd->pos1 + tsd_length - 1 < tir_pair->left_tir_start) {
          new_cost_left = tir_pair->left_tir_start
                            - (tsd->pos1 + tsd_length - 1);
        } else {
          new_cost_left = (tsd->pos1 + tsd_length - 1)
                            - tir_pair->left_tir_start;
        }

        if (tir_pair->right_tir_end < tsd->pos1 + tsd->offset) {
          new_cost_right = (tsd->pos1 + tsd->offset)
                             - tir_pair->right_tir_end;
        } else {
          new_cost_right = tir_pair->right_tir_end
                             - (tsd->pos1 + tsd->offset);
        }

        new_cost = new_cost_left + new_cost_right;
        if (new_cost < best_cost) {
          best_cost = new_cost;
          new_left_tir_start = tsd->pos1 + tsd_length;
          new_right_tir_end = tsd->pos1 + tsd->offset - 1;
          optimal_tsd_length = tsd_length;
        }
      }
    }
  }

  /* save the new borders and tsd length */
  if (info->TSDs.nextfreeSeed > 0) {
    tir_pair->left_tir_start = new_left_tir_start;
    tir_pair->right_tir_end = new_right_tir_end;
    tir_pair->tsd_length = optimal_tsd_length;
  }
}

/* this function searches for TSDs in the range of vicinity around the TIRs */
static int gt_tir_search_for_TSDs(GtTIRStream *tir_stream, TIRPair *tir_pair,
                                  const GtEncseq *encseq, GtError *err)
{
  unsigned long start_left_tir,
                end_left_tir,
                start_right_tir,
                end_right_tir,
                left_length,
                right_length,
                seq_end_pos,
                seq_start_pos,
                seq_length;
  unsigned long contignumber = tir_pair->contignumber;
  TSDinfo info;
  int had_err = 0;
  gt_error_check(err);

  /* check vicinity for left tir start */
  seq_start_pos = gt_encseq_seqstartpos(encseq, contignumber);
  seq_length = gt_encseq_seqlength(encseq, contignumber);

  gt_assert(tir_pair->left_tir_start >= seq_start_pos);
  gt_assert(tir_pair->left_tir_start <= tir_pair->left_tir_end);
  gt_assert(tir_pair->right_tir_start <= tir_pair->right_tir_end);
  gt_assert(tir_pair->right_tir_end <= seq_start_pos + seq_length);

  /* check if left tir start with vicinity aligns over sequence border */
  if (tir_pair->left_tir_start - seq_start_pos < tir_stream->vicinity)
    start_left_tir = seq_start_pos;
  else
    start_left_tir = tir_pair->left_tir_start - tir_stream->vicinity;

  /* do not align over end of left tir */
  if (tir_pair->left_tir_start + tir_stream->vicinity > tir_pair->left_tir_end)
    end_left_tir = tir_pair->left_tir_end;
  else
    end_left_tir = tir_pair->left_tir_start + tir_stream->vicinity;

  left_length = end_left_tir - start_left_tir + 1;

  /* vicinity of 3'-border of right tir
     do not align over 5'border of right tir */
  if (tir_pair->right_tir_end < tir_pair->right_tir_start
        + tir_stream->vicinity)
    start_right_tir = tir_pair->right_tir_start;
  else
    start_right_tir = tir_pair->right_tir_end - tir_stream->vicinity;

  seq_end_pos = seq_start_pos + seq_length - 1;

  /* do not align into next sequence in case of need decrease alignment
     length */
  if (tir_pair->right_tir_end + tir_stream->vicinity > seq_end_pos)
    end_right_tir = seq_end_pos;
  else
    end_right_tir = tir_pair->right_tir_end + tir_stream->vicinity;

  right_length = end_right_tir - start_right_tir + 1;

  /* search for TSDs */
  if (tir_stream->min_TSD_length > 1U) {
    /* dbseq (left) and query (right) are the encseqs wich will be aligned */
    GtUchar *dbseq, *query;
    dbseq = gt_calloc(left_length, sizeof (GtUchar));
    query = gt_calloc(right_length, sizeof (GtUchar));

    /* TODO: vor einem aufruf hier passiert ab einer gewissen
       größe ein speicherzugriffsfehler */
    gt_encseq_extract_decoded(encseq,(char*)dbseq,start_left_tir,end_left_tir);
    gt_encseq_extract_decoded(encseq,(char*)query,start_right_tir,
    end_right_tir);

    GT_INITARRAY(&info.TSDs, Seed);
    printf("TSD: %lu/%lu\n", start_left_tir, start_right_tir);
    printf("%s, %s\n", dbseq, query);
    gt_assert(start_left_tir < start_right_tir);
    info.left_start_pos = start_left_tir;
    info.right_start_pos = start_right_tir;

    if (gt_sarrquerysubstringmatch(dbseq,
                                   left_length,
                                   query,
                                   right_length,
                                   tir_stream->min_TSD_length,
                                   gt_encseq_alphabet(encseq),
                                   gt_tir_store_TSDs,
                                   &info,
                                   NULL,
                                   err) != 0) {
       had_err = -1;
    }
    gt_free(dbseq);
    gt_free(query);

    /* find the best TSD */
    if (!had_err)
      gt_tir_find_best_TSD(&info, tir_stream, tir_pair);

    GT_FREEARRAY(&info.TSDs, Seed);
  }

  return had_err;
}

static int gt_tir_searchforTIRs(GtTIRStream *tir_stream,
                                GtArrayTIRPair *tir_pairs,
                                const GtEncseq *encseq, GtError *err)
{
  int count;
  unsigned long seedcounter = 0;
  GtArrayTIRPair new;             /* need to remove overlaps */
  unsigned long right_tir_start;  /* need to calculate the reverse position */
  unsigned long right_tir_end;    /* need to calculate the reverse position */
  GtXdropresources *xdropresources;
  unsigned long total_length = 0,
                left_tir_length = 0,
                right_tir_length = 0,
                max_left_length = 0,    /* maximal length of left TIR*/
                max_right_length = 0,
                distance = 0;
  GtUchar *left_tir_char = NULL,        /* next character to align */
          *right_tir_char = NULL;
  unsigned long edist = 0,
                alilen,
                seqstart, seqend;
  Seed *seedptr;
  TIRPair *pair;
  int had_err = 0;
  GtXdropbest xdropbest_left, xdropbest_right;
  GtSeqabstract *sa_useq = gt_seqabstract_new_empty(),
                *sa_vseq = gt_seqabstract_new_empty();
  gt_error_check(err);

  printf("searching for TIRs\n");

  xdropresources = gt_xdrop_resources_new(&tir_stream->arbit_scores);

  /* Iterating over seeds */
  for (seedcounter = 0; seedcounter < tir_stream->seedinfo.seed.nextfreeSeed;
       seedcounter++) {
    GtUchar *dbseq, *query;

    seedptr = &(tir_stream->seedinfo.seed.spaceSeed[seedcounter]);
    gt_assert(tir_stream->seedinfo.max_tir_length >= seedptr->len);
    alilen = tir_stream->seedinfo.max_tir_length - seedptr->len;
    seqstart = gt_encseq_seqstartpos(tir_stream->encseq, seedptr->contignumber);
    seqend = seqstart
                + gt_encseq_seqlength(tir_stream->encseq,seedptr->contignumber);

    /* left (reverse) xdrop */
    if (seedptr->pos1 > 0)
    {
      if (alilen <= seedptr->pos1)
      {
        gt_seqabstract_reinit_encseq(sa_useq, encseq, alilen, 0);
        gt_seqabstract_reinit_encseq(sa_vseq, encseq, alilen, 0);
      } else
      {
        unsigned long maxleft = seedptr->pos1 - seqstart;
        gt_seqabstract_reinit_encseq(sa_useq, encseq, maxleft, 0);
        gt_seqabstract_reinit_encseq(sa_vseq, encseq, maxleft, 0);
      }
      gt_evalxdroparbitscoresextend(false,
                                    &xdropbest_left,
                                    xdropresources,
                                    sa_useq,
                                    sa_vseq,
                                    seedptr->pos1,
                                    seedptr->pos2 + seedptr->offset,
                                   (GtXdropscore) tir_stream->xdrop_belowscore);
    } else
    {
      xdropbest_left.ivalue = 0;
      xdropbest_left.jvalue = 0;
      xdropbest_left.score = 0;
    }
/*
    gt_evalxdroparbitscoresright(&tir_stream->arbit_scores,
                                &xdropbest_right,
                                &fronts,
                                encseq,
                                encseq,
                                seedptr->pos1 + seedptr->len,
                                seedptr->pos2 + seedptr->len,
                                (int) (total_length -
                                                (seedptr->pos1 + seedptr->len)),
                                (int) (total_length -
                                                (seedptr->pos2 + seedptr->len)),
                                (Xdropscore)tir_stream->xdrop_belowscore);
    GT_FREEARRAY (&fronts, Myfrontvalue); */

    total_length = gt_encseq_total_length(encseq);
    if (seedptr->pos1 + seedptr->len > 0)
    {
      gt_assert(seqend >= seedptr->pos1 + seedptr->offset + seedptr->len);
      if (alilen <= seqend - (seedptr->pos1 + seedptr->offset +
                              seedptr->len))
      {
        gt_seqabstract_reinit_encseq(sa_useq,tir_stream->encseq,alilen,0);
        gt_seqabstract_reinit_encseq(sa_vseq,tir_stream->encseq,alilen,0);
      } else
      {
        gt_seqabstract_reinit_encseq(sa_useq,tir_stream->encseq,
                                     seqend - (seedptr->pos1+seedptr->len),
                                     0);
        gt_seqabstract_reinit_encseq(sa_vseq,tir_stream->encseq,
                                     seqend - (seedptr->pos1 +
                                               seedptr->offset +
                                               seedptr->len),0);
      }
      gt_evalxdroparbitscoresextend(true,
                                    &xdropbest_right,
                                    xdropresources,
                                    sa_useq,
                                    sa_vseq,
                                    seedptr->pos1 + seedptr->len,
                                    seedptr->pos1 + seedptr->offset +
                                    seedptr->len,
                                   (GtXdropscore) tir_stream->xdrop_belowscore);
    } else
    {
      xdropbest_right.ivalue = 0;
      xdropbest_right.jvalue = 0;
      xdropbest_right.score = 0;
    }

    GT_GETNEXTFREEINARRAY(pair, tir_pairs, TIRPair, 5);
    /* Store positions for the found TIR */
    pair->contignumber = seedptr->contignumber;
    pair->left_tir_start = seedptr->pos1 - xdropbest_left.ivalue;
    pair->left_tir_end = seedptr->pos1 + seedptr->len - 1
                            + xdropbest_right.ivalue;
printf("%lu %lu %lu\n", pair->contignumber, pair->left_tir_start,
                     pair->left_tir_end);

    /* We want the actual positions (not mirrored)
       end of mirrored is start of actual TIR */
    right_tir_start = seedptr->pos2 + seedptr->len - 1 + xdropbest_right.jvalue;
    /* start of mirrored is end of actual TIR */
    right_tir_end   = seedptr->pos2 - xdropbest_left.jvalue;

    /* We can get the corresponding position of a mirrored one with
       GT_REVERSEPOS(total length,position) */
    pair->right_tir_start = GT_REVERSEPOS(total_length,right_tir_start);
    pair->right_tir_end = GT_REVERSEPOS(total_length,right_tir_end);
    pair->similarity = 0.0;
    pair->skip = false;
  }

  /* initialize array for removing overlaps */
  GT_INITARRAY(&new, TIRPair);

  /* remove overlaps if wanted */
  if (tir_stream->best_overlaps || tir_stream->no_overlaps) {
    tir_stream->num_of_tirs = gt_tir_remove_overlaps(tir_pairs, &new,
                                                     tir_stream->num_of_tirs,
                                                     tir_stream->no_overlaps);

    /* set references to the new array and free old array */
    GT_FREEARRAY(tir_pairs, TIRPair);
    tir_stream->tir_pairs = new;
    tir_pairs = &new;
  }

  /* sort the tir_pairs */
  gt_tir_sort_TIRs(tir_pairs, 0, tir_stream->num_of_tirs);

  /* log output on console */
  gt_log_log("TIRs found:\n");
  for (count = 0; count < tir_stream->tir_pairs.nextfreeTIRPair; count++) {
    gt_log_log("contig %lu\n left: start %lu, end %lu\n right: start %lu, "
               "end %lu\n similarity: %f\n\n",
               tir_stream->tir_pairs.spaceTIRPair[count].contignumber,
               tir_stream->tir_pairs.spaceTIRPair[count].left_tir_start,
               tir_stream->tir_pairs.spaceTIRPair[count].left_tir_end,
               tir_stream->tir_pairs.spaceTIRPair[count].right_tir_start,
               tir_stream->tir_pairs.spaceTIRPair[count].right_tir_end,
               tir_stream->tir_pairs.spaceTIRPair[count].similarity);
  }

  return had_err;
}

static int gt_tir_stream_next(GtNodeStream *ns, GT_UNUSED GtGenomeNode **gn,
                              GtError *err)
{
  GtTIRStream *tir_stream;
  int had_err = 0;
  gt_error_check(err);
  tir_stream = gt_node_stream_cast(gt_tir_stream_class(), ns);

  /* generate and check seeds */
   if (tir_stream->state == GT_TIR_STREAM_STATE_START) {
    if (!had_err && gt_enumeratemaxpairs(tir_stream->ssar,
                      tir_stream->encseq,
                      gt_readmodeSequentialsuffixarrayreader(tir_stream->ssar),
                      (unsigned int) tir_stream->min_seed_length,
                      gt_tir_store_seeds,
                      &tir_stream->seedinfo,
                      err) != 0) {
      had_err = -1;
    }

    /* extend seeds to TIRs and check TIRs */
    if (!had_err && gt_tir_searchforTIRs(tir_stream,
                                         &tir_stream->tir_pairs,
                                         tir_stream->encseq, err) != 0) {
      had_err = -1;
    }

    GT_FREEARRAY(&tir_stream->seedinfo.seed, Seed);
    tir_stream->state = GT_TIR_STREAM_STATE_REGIONS;
  }

  /* stream out the region nodes */
  if (!had_err && tir_stream->state == GT_TIR_STREAM_STATE_REGIONS) {
    bool skip = false;

    /* check whether index is valid */
    if (tir_stream->cur_elem_index < tir_stream->num_of_tirs) {
      unsigned long seqnum,
                    seqlength;
      GtGenomeNode *rn;
      GtStr *seqid;
      seqnum = tir_stream->tir_pairs.spaceTIRPair[tir_stream->cur_elem_index]
                                    .contignumber;

      /* if first time we do this */
      if (tir_stream->prev_seqnum == GT_UNDEF_ULONG) {
        /* use current seqnum */
        tir_stream->prev_seqnum = seqnum;
      } else {
        /* else get seqnum of next contig */
        while (tir_stream->prev_seqnum == seqnum) {
          tir_stream->cur_elem_index++;

          /* do not go on if index is out of bounds */
          if (tir_stream->cur_elem_index >= tir_stream->num_of_tirs) {
            skip = true;
            break;
          }

          seqnum = tir_stream->tir_pairs
                                       .spaceTIRPair[tir_stream->cur_elem_index]
                                       .contignumber;
        }
      }

      /* create node */
      if (!skip) {
        tir_stream->prev_seqnum = seqnum;
        seqlength = gt_encseq_seqlength(tir_stream->encseq, seqnum);

        /* TODO: add real description */
        seqid = gt_str_new_cstr("seq");
        gt_str_append_ulong(seqid, seqnum);
        rn = gt_region_node_new(seqid, 1, seqlength);

        gt_str_delete(seqid);
        *gn = rn;
        tir_stream->cur_elem_index++;
      } else {
        /* skipping */
        tir_stream->cur_elem_index = 0;
        tir_stream->state = GT_TIR_STREAM_STATE_COMMENTS;
        *gn = NULL;
      }
    } else {
      /* no valid index */
      tir_stream->cur_elem_index = 0;
      tir_stream->state = GT_TIR_STREAM_STATE_COMMENTS;
      *gn = NULL;
    }
  }

  /* then stream out the comment nodes */
  if (!had_err && tir_stream->state == GT_TIR_STREAM_STATE_COMMENTS)
  {
    bool skip = false;
    if (tir_stream->cur_elem_index < tir_stream->num_of_tirs) {
      const char *description;
      unsigned long description_len, seqnum;
      GtGenomeNode *cn;

      seqnum = tir_stream->tir_pairs.spaceTIRPair[tir_stream->cur_elem_index]
                                    .contignumber;

      /* for the first time */
      if (tir_stream->prev_seqnum == GT_UNDEF_LONG) {
        /* use current seqnum */
        tir_stream->prev_seqnum = seqnum;
      } else {
        /* else get seqnum of next contig */
        while (tir_stream->prev_seqnum == seqnum) {
          tir_stream->cur_elem_index++;
          if (tir_stream->cur_elem_index >= tir_stream->num_of_tirs) {
            skip = true;
            break;
          }
        }

        seqnum = tir_stream->tir_pairs.spaceTIRPair[tir_stream->cur_elem_index]
                                      .contignumber;
      }

      /* create a new comment node */
      if (!skip) {
        char description_string[BUFSIZ];
        tir_stream->prev_seqnum = seqnum;

        /* get description and descriptionlength of current contig */
        description = gt_encseq_description(tir_stream->encseq,
                                            &description_len, seqnum);

       (void) strncpy(description_string, description,
                      (size_t) (description_len * sizeof (char)));
        description_string[description_len] = '\0';

        /* make a new comment node */
        cn = gt_comment_node_new(description_string);

        *gn = cn;
        tir_stream->cur_elem_index++;
      } else {
        /* skipping */
        tir_stream->cur_elem_index = 0;
        tir_stream->state = GT_TIR_STREAM_STATE_FEATURES;
        *gn = NULL;
      }
    } else {
      /* no valid index */
      tir_stream->cur_elem_index = 0;
      tir_stream->state = GT_TIR_STREAM_STATE_FEATURES;
      *gn = NULL;
    }
  }

  /* finally stream out the features */
  if (!had_err && tir_stream->state == GT_TIR_STREAM_STATE_FEATURES) {
    if (tir_stream->cur_elem_index < tir_stream->num_of_tirs) {
      GtStr *seqid, *source;
      GtGenomeNode *node, GT_UNUSED *parent;
      const TIRPair *pair =
                &tir_stream->tir_pairs.spaceTIRPair[tir_stream->cur_elem_index];
      unsigned long seqstartpos;

      seqstartpos = gt_encseq_seqstartpos(tir_stream->encseq,
                                          pair->contignumber);
      seqid = gt_str_new_cstr("seq");
      source = gt_str_new_cstr("TIRvish");

      gt_str_append_ulong(seqid, pair->contignumber);

      /* repeat region */

      node = gt_feature_node_new(seqid, "repeat_region",
                                 pair->left_tir_start - seqstartpos + 1,
                                 pair->right_tir_end - seqstartpos + 1,
                                 GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*) node, source);
      *gn = node;
      parent = node;

      /* target site duplication */

      if (tir_stream->min_TSD_length > 1U) {
        node = gt_feature_node_new(seqid,
                       "target_site_duplication",
                       pair->left_tir_start - seqstartpos + 1
                         - pair->tsd_length,
                       pair->left_tir_start - seqstartpos,
                       GT_STRAND_UNKNOWN);
        gt_feature_node_set_source((GtFeatureNode*) node, source);
        gt_feature_node_add_child((GtFeatureNode*) parent,
                      (GtFeatureNode*) node);
        node = gt_feature_node_new(seqid,
                       "target_site_duplication",
                       pair->right_tir_end - seqstartpos + 1, /* XXX: why was
                                                                 here +2 in
                                                                 LTRharvest? */
                       pair->right_tir_end - seqstartpos + 1
                       + pair->tsd_length,
                       GT_STRAND_UNKNOWN);
        gt_feature_node_set_source((GtFeatureNode*) node, source);
        gt_feature_node_add_child((GtFeatureNode*) parent,
                      (GtFeatureNode*) node);
      }

      /* terminal inverted repeat element */

      node = gt_feature_node_new(seqid, "terminal_inverted_repeat_element",
                                 pair->left_tir_start - seqstartpos + 1,
                                 pair->right_tir_end - seqstartpos +1,
                                 GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*)node, source);
      gt_feature_node_add_child((GtFeatureNode*)parent, (GtFeatureNode*)node);
      parent = node;

      /* left terminal inverted repeat */

      node = gt_feature_node_new(seqid, "terminal_inverted_repeat",
                                 pair->left_tir_start - seqstartpos + 1,
                                 pair->left_tir_end - seqstartpos + 1,
                                 GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*)node, source);
      gt_feature_node_add_child((GtFeatureNode*)parent, (GtFeatureNode*)node);

      /* right terminal inverted repeat */

      node = gt_feature_node_new(seqid, "terminal_inverted_repeat",
                                 pair->right_tir_start - seqstartpos + 1,
                                 pair->right_tir_end - seqstartpos + 1,
                                 GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*)node, source);
      gt_feature_node_add_child((GtFeatureNode*)parent, (GtFeatureNode*)node);

      /* clean up and get next pair */
      gt_str_delete(seqid);
      gt_str_delete(source);
      tir_stream->cur_elem_index++;
    } else {
      *gn = NULL;
    }
  }
  return had_err;
}

static void gt_tir_stream_free(GtNodeStream *ns)
{
  GtTIRStream *tir_stream = gt_node_stream_cast(gt_tir_stream_class(), ns);
  gt_str_delete(tir_stream->str_indexname);
  if (tir_stream->ssar != NULL)
    gt_freeSequentialsuffixarrayreader(&tir_stream->ssar);
}

const GtNodeStreamClass* gt_tir_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtTIRStream),
                                   gt_tir_stream_free,
                                   gt_tir_stream_next);
  }
  return nsc;
}

GtNodeStream* gt_tir_stream_new(GtStr *str_indexname,
                                unsigned long min_seed_length,
                                unsigned long min_TIR_length,
                                unsigned long max_TIR_length,
                                unsigned long min_TIR_distance,
                                unsigned long max_TIR_distance,
                                GtXdropArbitraryscores arbit_scores,
                                int xdrop_belowscore,
                                double similarity_threshold,
                                bool best_overlaps,
                                bool no_overlaps,
                                unsigned long min_TSD_length,
                                unsigned long max_TSD_length,
                                unsigned long vicinity,
                                GtError *err)
{
  int had_err = 0;
  GtNodeStream *ns = gt_node_stream_create(gt_tir_stream_class(), false);
  GtTIRStream *tir_stream = gt_node_stream_cast(gt_tir_stream_class(), ns);
  tir_stream->num_of_tirs = 0;
  tir_stream->cur_elem_index = 0;
  tir_stream->prev_seqnum = GT_UNDEF_ULONG;
  tir_stream->state = GT_TIR_STREAM_STATE_START;

  tir_stream->str_indexname = gt_str_ref(str_indexname);
  tir_stream->min_seed_length = min_seed_length;
  tir_stream->min_TIR_length = min_TIR_length;
  tir_stream->max_TIR_length = max_TIR_length;
  tir_stream->min_TIR_distance = min_TIR_distance;
  tir_stream->max_TIR_distance = max_TIR_distance;
  tir_stream->arbit_scores = arbit_scores;
  tir_stream->xdrop_belowscore = xdrop_belowscore;
  tir_stream->similarity_threshold = similarity_threshold;
  tir_stream->best_overlaps = best_overlaps;
  tir_stream->no_overlaps = no_overlaps;
  tir_stream->min_TSD_length = min_TSD_length;
  tir_stream->max_TSD_length = max_TSD_length;
  tir_stream->vicinity = vicinity;

  tir_stream->seedinfo.max_tir_length = max_TIR_length;
  tir_stream->seedinfo.min_tir_length = min_TIR_length;
  tir_stream->seedinfo.max_tir_distance = max_TIR_distance;
  tir_stream->seedinfo.min_tir_distance = min_TIR_distance;

  tir_stream->ssar =
      gt_newSequentialsuffixarrayreaderfromfile(gt_str_get(str_indexname),
                                                SARR_LCPTAB | SARR_SUFTAB |
                                                SARR_ESQTAB | SARR_DESTAB |
                                                SARR_SSPTAB | SARR_SDSTAB,
                                                SEQ_scan,
                                                NULL,
                                                err);
  if (tir_stream->ssar == NULL) {
    gt_node_stream_delete(ns);
    return NULL;
  }
  tir_stream->encseq = gt_encseqSequentialsuffixarrayreader(tir_stream->ssar);
  gt_assert(gt_encseq_is_mirrored(tir_stream->encseq));
  GT_INITARRAY(&tir_stream->seedinfo.seed, Seed);
  GT_INITARRAY(&tir_stream->tir_pairs, TIRPair);
  tir_stream->seedinfo.num_of_contigs =
                                 gt_encseq_num_of_sequences(tir_stream->encseq);
  tir_stream->seedinfo.totallength =
                                 gt_encseq_total_length(tir_stream->encseq);
  tir_stream->seedinfo.midpos =
                           (gt_encseq_total_length(tir_stream->encseq) - 1) / 2;

  if (!had_err)
    return ns;

  return NULL;
}
