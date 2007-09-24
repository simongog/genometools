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

#include <inttypes.h>
#include "libgtcore/strarray.h"
#include "symboldef.h"
#include "arraydef.h"
#include "chardef.h"
#include "fbs-def.h"
#include "seqdesc.h"
#include "format64.h"
#include "iterseq.h"

#include "fbsadv.pr"

#include "readnextUchar.gen"

int overallquerysequences(int(*processsequence)(void *,
                                                uint64_t,
                                                const Uchar *,
                                                unsigned long,
                                                const char *,
                                                Env *),
                          void *info,
                          ArrayUchar *sequencebuffer,
                          const StrArray *filenametab,
                          Sequencedescription *sequencedescription,
                          const Uchar *symbolmap,
                          Env *env)
{
  Fastabufferstate fbs;
  Uchar charcode;
  int retval;
  uint64_t unitnum = 0;
  bool haserr = false;
  char *desc;

  env_error_check(env);
  initformatbufferstate(&fbs,
                        filenametab,
                        symbolmap,
                        false,
                        NULL,
                        sequencedescription,
                        NULL,
                        env);
  sequencebuffer->nextfreeUchar = 0;
  while (true)
  {
    retval = readnextUchar(&charcode,&fbs,env);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      break;
    }
    if (charcode == (Uchar) SEPARATOR)
    {
      if (sequencebuffer->nextfreeUchar == 0)
      {
        env_error_set(env,"sequence " Formatuint64_t " is empty",
                      PRINTuint64_tcast(unitnum));
        haserr = true;
        break;
      }
      desc = queue_get(sequencedescription->descptr,env);
      if (processsequence(info,
                          unitnum,
                          sequencebuffer->spaceUchar,
                          sequencebuffer->nextfreeUchar,
                          desc,
                          env) != 0)
      {
        haserr = true;
        FREESPACE(desc);
        break;
      }
      FREESPACE(desc);
      sequencebuffer->nextfreeUchar = 0;
      unitnum++;
    } else
    {
      STOREINARRAY(sequencebuffer,Uchar,1024,charcode);
    }
  }
  if (!haserr && sequencebuffer->nextfreeUchar > 0)
  {
    desc = queue_get(sequencedescription->descptr,env);
    if (processsequence(info,
                        unitnum,
                        sequencebuffer->spaceUchar,
                        sequencebuffer->nextfreeUchar,
                        desc,
                        env) != 0)
    {
      haserr = true;
    }
    FREESPACE(desc);
    sequencebuffer->nextfreeUchar = 0;
  }
  return haserr ? -1 : 0;
}

 struct Scansequenceiterator
{
  Fastabufferstate fbs;
  const StrArray *filenametab;
  const Uchar *symbolmap;
  Sequencedescription sequencedescription;
  ArrayUchar sequencebuffer;
  uint64_t unitnum;
  bool withsequence, exhausted;
};

Scansequenceiterator *newScansequenceiterator(const StrArray *filenametab,
                                              const Uchar *symbolmap,
                                              bool withsequence,
                                              Env *env)
{
  Scansequenceiterator *sseqit;

  ALLOCASSIGNSPACE(sseqit,NULL,Scansequenceiterator,1);
  INITARRAY(&sseqit->sequencebuffer,Uchar);
  INITARRAY(&sseqit->sequencedescription.headerbuffer,char);
  sseqit->sequencedescription.descptr = queue_new(env);
  initformatbufferstate(&sseqit->fbs,
                        filenametab,
                        symbolmap,
                        false,
                        NULL,
                        &sseqit->sequencedescription,
                        NULL,
                        env);
  sseqit->sequencebuffer.nextfreeUchar = 0;
  sseqit->exhausted = false;
  sseqit->unitnum = 0;
  sseqit->withsequence = withsequence;
  return sseqit;
}

void freeScansequenceiterator(Scansequenceiterator **sseqit,Env *env)
{
  queue_delete_with_contents((*sseqit)->sequencedescription.descptr,env);
  FREEARRAY(&(*sseqit)->sequencebuffer,Uchar);
  FREEARRAY(&(*sseqit)->sequencedescription.headerbuffer,char);
  FREESPACE(*sseqit);
}

int nextScansequenceiterator(const Uchar **sequence,
                             unsigned long *len,
                             char **desc,
                             Scansequenceiterator *sseqit,
                             Env *env)
{
  Uchar charcode;
  int retval;
  bool haserr = false, foundseq = false;

  if (sseqit->exhausted)
  {
    return 0;
  }
  while (true)
  {
    retval = readnextUchar(&charcode,&sseqit->fbs,env);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      sseqit->exhausted = true;
      break;
    }
    if (charcode == (Uchar) SEPARATOR)
    {
      if (sseqit->sequencebuffer.nextfreeUchar == 0 && sseqit->withsequence)
      {
        env_error_set(env,"sequence " Formatuint64_t " is empty",
                      PRINTuint64_tcast(sseqit->unitnum));
        haserr = true;
        break;
      }
      *desc = queue_get(sseqit->sequencedescription.descptr,env);
      *len = sseqit->sequencebuffer.nextfreeUchar;
      if (sseqit->withsequence)
      {
        *sequence = sseqit->sequencebuffer.spaceUchar;
      }
      sseqit->sequencebuffer.nextfreeUchar = 0;
      foundseq = true;
      sseqit->unitnum++;
      break;
    } else
    {
      if (sseqit->withsequence)
      {
        STOREINARRAY(&sseqit->sequencebuffer,Uchar,1024,charcode);
      } else
      {
        sseqit->sequencebuffer.nextfreeUchar++;
      }
    }
  }
  if (!haserr && sseqit->sequencebuffer.nextfreeUchar > 0)
  {
    *desc = queue_get(sseqit->sequencedescription.descptr,env);
    if (sseqit->withsequence)
    {
      *sequence = sseqit->sequencebuffer.spaceUchar;
    }
    *len = sseqit->sequencebuffer.nextfreeUchar;
    foundseq = true;
    sseqit->sequencebuffer.nextfreeUchar = 0;
  }
  if (haserr)
  {
    return -1;
  }
  if (foundseq)
  {
    return 1;
  }
  return 0;
}
