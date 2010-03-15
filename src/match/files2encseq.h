/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef FILES2ENCSEQ_H
#define FILES2ENCSEQ_H

#include <stdbool.h>
#include "core/str_api.h"
#include "core/str_array_api.h"
#include "seqpos-def.h"
#include "sfx-progress.h"
#include "encseq-def.h"
#include "core/logger.h"

Encodedsequence *fromfiles2encseq(ArraySeqpos *sequenceseppos,
                                  Sfxprogress *sfxprogress,
                                  const GtStr *str_indexname,
                                  const GtStr *str_smap,
                                  const GtStr *str_sat,
                                  const GtStrArray *filenametab,
                                  bool isdna,
                                  bool isprotein,
                                  bool isplain,
                                  bool outtistab,
                                  bool outdestab,
                                  bool outsdstab,
                                  bool outssptab,
                                  GtLogger *logger,
                                  GtError *err);

#endif
