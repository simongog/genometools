/*
  Copyright (c) 2007-2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Center for Bioinformatics, University of Hamburg

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

#ifndef ESA_MAXPAIRS_H
#define ESA_MAXPAIRS_H

#include "core/error_api.h"
#include "seqpos-def.h"
#include "encseq-def.h"
#include "esa-seqread.h"
#include "core/logger.h"

typedef int (*Processmaxpairs)(void *,const Encodedsequence *,
                               Seqpos,Seqpos,Seqpos,GtError *);

int enumeratemaxpairs(Sequentialsuffixarrayreader *ssar,
                      const Encodedsequence *encseq,
                      Readmode readmode,
                      unsigned int searchlength,
                      Processmaxpairs processmaxpairs,
                      void *processmaxpairsinfo,
                      GtLogger *logger,
                      GtError *err);

#endif
