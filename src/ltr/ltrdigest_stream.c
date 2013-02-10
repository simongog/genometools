/*
  Copyright (c) 2008-2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2009 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/class_alloc_lock.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/range.h"
#include "core/str.h"
#include "core/symbol.h"
#include "core/unused_api.h"
#include "extended/extract_feature_sequence.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_stream_api.h"
#include "extended/reverse_api.h"
#include "ltr/ltrdigest_def.h"
#include "ltr/ltrdigest_stream.h"
#include "ltr/ltr_visitor.h"

struct GtLTRdigestStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtEncseq *encseq;
  GtRegionMapping *rmap;
  GtLTRVisitor *lv;
  GtStr *ltrdigest_tag;
  int tests_to_run;
  GtLTRElement element;
};

#define GT_ALIWIDTH 60

#define gt_ltrdigest_stream_cast(GS)\
        gt_node_stream_cast(gt_ltrdigest_stream_class(), GS)

static int run_ltrdigest(GtLTRElement *element, char *seq,
                         GT_UNUSED GtLTRdigestStream *ls,
#ifdef HAVE_HMMER
                         GtError *err)
#else
                         GT_UNUSED GtError *err)
#endif
{
  int had_err = 0;
  char *rev_seq;
  unsigned long seqlen = gt_ltrelement_length(element);
  GT_UNUSED GtStrand canonical_strand = GT_STRAND_UNKNOWN;

  /* create reverse strand sequence */
  rev_seq = gt_calloc((size_t) seqlen+1, sizeof (char));
  memcpy(rev_seq, seq, sizeof (char) * seqlen);
  had_err = gt_reverse_complement(rev_seq, seqlen, err);

  if (!had_err)
  {
#ifdef HAVE_HMMER
    /* Protein domain finding
       ---------------------- */
    /*if (ls->tests_to_run & GT_LTRDIGEST_RUN_PDOM)
    {
      GtPdomResults *pdom_results = NULL;
      if (!ls->pdf)
      {
        gt_error_set(err, "No PdomFinder object found -- "
                          "how could that happen?");
        had_err = -1;
      } else
      {
        pdom_results = gt_pdom_finder_find(ls->pdf, (const char*) seq,
                                           (const char*) rev_seq, element, err);
        if (!pdom_results)
        {
          had_err = -1;
        } else {
          if (pdom_results && !gt_pdom_results_empty(pdom_results))
          {
            if (gt_double_compare(
                     gt_pdom_results_get_combined_evalue_fwd(pdom_results),
                     gt_pdom_results_get_combined_evalue_rev(pdom_results)) < 0)
              canonical_strand = GT_STRAND_FORWARD;
            else
              canonical_strand = GT_STRAND_REVERSE;
            gt_feature_node_set_strand(ls->element.mainnode, canonical_strand);
            (void) gt_pdom_results_foreach_domain_hit(pdom_results,
                                                      pdom_hit_attach_gff3,
                                                      ls,
                                                      err);
          }
          gt_pdom_results_delete(pdom_results);
        }
      }
    }*/
#endif

    /* PPT finding
       ----------- */
    /* if (ls->tests_to_run & GT_LTRDIGEST_RUN_PPT)
    {
      GtPPTResults *ppt_results = NULL;
      ppt_results = gt_ppt_find((const char*) seq, (const char*) rev_seq,
                              element, ls->ppt_opts);
      if (gt_ppt_results_get_number_of_hits(ppt_results) > 0)
      {
        ppt_attach_results_to_gff3(ppt_results, element, &canonical_strand,
                                   ls->ltrdigest_tag);
      }
      gt_ppt_results_delete(ppt_results);
    }*/

    /* PBS finding
       ----------- */
    /* if (ls->tests_to_run & GT_LTRDIGEST_RUN_PBS)
    {
      GtPBSResults *pbs_results = NULL;
      pbs_results = gt_pbs_find((const char*) seq, (const char*) rev_seq,
                                element, ls->pbs_opts, err);
       if (gt_pbs_results_get_number_of_hits(pbs_results) > 0)
       {
        pbs_attach_results_to_gff3(pbs_results, element, &canonical_strand,
                                   ls->ltrdigest_tag);
       }
       gt_pbs_results_delete(pbs_results);
    }*/
  }
  gt_free(rev_seq);
  return had_err;
}

static void set_element_direction(GtFeatureNode *rootnode, GtStrand strand)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *fn = rootnode;
  fni = gt_feature_node_iterator_new(rootnode);
  for (fn = rootnode; fn != NULL; fn = gt_feature_node_iterator_next(fni)) {
    GtStrand cur_strand = gt_feature_node_get_strand(fn);
    if (cur_strand == GT_STRAND_UNKNOWN) {
      gt_feature_node_set_strand(fn, strand);
    }
  }
  gt_feature_node_iterator_delete(fni);
}

static int gt_ltrdigest_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                    GtError *err)
{
  GtLTRdigestStream *ls;
  GtFeatureNode *fn = NULL;
  int had_err;

  gt_error_check(err);
  ls = gt_ltrdigest_stream_cast(ns);

  /* initialize this element */
  memset(&ls->element, 0, sizeof (GtLTRElement));

  /* get annotations from parser */
  had_err = gt_node_stream_next(ls->in_stream, gn, err);
  if (!had_err && *gn)
  {
    GtFeatureNodeIterator *gni;
    GtFeatureNode *mygn;

   /* only process feature nodes */
   if (!(fn = gt_feature_node_try_cast(*gn)))
     return 0;

    /* fill LTRElement structure from GFF3 subgraph */
    gni = gt_feature_node_iterator_new(fn);
    for (mygn = fn; mygn != NULL; mygn = gt_feature_node_iterator_next(gni))
      (void) gt_genome_node_accept((GtGenomeNode*) mygn,
                                   (GtNodeVisitor*) ls->lv,
                                   err);
    gt_feature_node_iterator_delete(gni);
  }

  if (ls->element.mainnode != NULL)
  {
    GtStr *seq;

    seq = gt_str_new();
    had_err = gt_extract_feature_sequence(seq,
                                          (GtGenomeNode*) ls->element.mainnode,
                                          gt_symbol(gt_ft_LTR_retrotransposon),
                                          false, NULL, NULL, ls->rmap, err);

    if (!had_err) {
      /* run LTRdigest core routine */
      had_err = run_ltrdigest(&ls->element, gt_str_get(seq), ls, err);
    }
    gt_str_delete(seq);
    if (!had_err) {
      GtStrand ltr_strand = gt_feature_node_get_strand(ls->element.mainnode);
      if (ltr_strand != GT_STRAND_UNKNOWN) {
        set_element_direction(fn, ltr_strand);
      }
    }
  }
  if (had_err) {
    gt_genome_node_delete(*gn);
    *gn = NULL;
  }
  return had_err;
}

static void gt_ltrdigest_stream_free(GtNodeStream *ns)
{
  GtLTRdigestStream *ls = gt_ltrdigest_stream_cast(ns);
  gt_node_visitor_delete((GtNodeVisitor*) ls->lv);
  gt_str_delete(ls->ltrdigest_tag);
  gt_node_stream_delete(ls->in_stream);
#ifdef HAVE_HMMER
  /* gt_pdom_finder_delete(ls->pdf); */
#endif
}

const GtNodeStreamClass* gt_ltrdigest_stream_class(void)
{
  static const GtNodeStreamClass *nsc;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtLTRdigestStream),
                                   gt_ltrdigest_stream_free,
                                   gt_ltrdigest_stream_next );
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_ltrdigest_stream_new(GtNodeStream *in_stream,
                                      int tests_to_run,
                                      GtRegionMapping *rmap,
                                      GT_UNUSED GtError *err)
{
  GtNodeStream *ns;
  GtLTRdigestStream *ls;
  ns = gt_node_stream_create(gt_ltrdigest_stream_class(), true);
  ls = gt_ltrdigest_stream_cast(ns);
  ls->in_stream = gt_node_stream_ref(in_stream);
  ls->tests_to_run = tests_to_run;
  ls->rmap = rmap;
  ls->ltrdigest_tag = gt_str_new_cstr(GT_LTRDIGEST_TAG);
  ls->lv = (GtLTRVisitor*) gt_ltr_visitor_new(&ls->element);
  return ns;
}
