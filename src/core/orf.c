/*
  Copyright (c) 2006, 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006       Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include "core/codon.h"
#include "core/orf.h"
#include "core/range.h"
#include "core/undef.h"

#define START_AMINO  'M'
#define STOP_AMINO   '*'

void gt_determine_ORFs(GtORFProcessor orf_processor, void *data,
                       unsigned int framenum, const char *frame,
                       unsigned long framelen, bool start_codon,
                       bool final_stop_codon)
{
  unsigned long i;
  GtRange orf;
  gt_assert(orf_processor && framenum <= 2 && frame);
  orf.start = GT_UNDEF_ULONG;
  for (i = 0; i < framelen; i++) {
    if (orf.start == GT_UNDEF_ULONG) {
      if (start_codon) {
        if (frame[i] == START_AMINO)
          orf.start = i * GT_CODON_LENGTH + framenum;
      }
      else {
        if (i) {
          if (frame[i-1] == STOP_AMINO && frame[i] != STOP_AMINO)
            orf.start = i * GT_CODON_LENGTH + framenum;
        }
        else {
          if (frame[i] != STOP_AMINO)
            orf.start = framenum;
        }
      }
    }
    else {
      if (frame[i] == STOP_AMINO) {
        orf.end = i * GT_CODON_LENGTH + framenum + 2;
        orf_processor(data, &orf);
        orf.start = GT_UNDEF_ULONG;
      }
    }
  }
  if (!final_stop_codon) {
    if (orf.start != GT_UNDEF_ULONG) {
      orf.end = (framelen - 1) * GT_CODON_LENGTH + framenum + 2;
      orf_processor(data, &orf);
    }
  }
}
