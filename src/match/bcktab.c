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

#include <string.h>
#include <errno.h>
#include <math.h>
#include "core/error.h"
#include "core/str.h"
#include "core/fa.h"
#include "core/types_api.h"
#include "core/chardef.h"
#include "core/mathsupport.h"
#include "core/unused_api.h"
#include "core/format64.h"
#include "core/mapspec-gen.h"
#include "core/unused_api.h"
#include "core/minmax.h"
#include "core/log_api.h"
#include "core/safecast-gen.h"
#include "intcode-def.h"
#include "esa-fileend.h"
#include "bcktab.h"
#include "initbasepower.h"
#include "stamp.h"

#define FROMCODE2SPECIALCODE(CODE,NUMOFCHARS)\
                            (((CODE) - ((NUMOFCHARS)-1)) / (NUMOFCHARS))

typedef struct
{
  unsigned long nonspecialsmaxbucketsize,
                specialsmaxbucketsize,
                maxbucketsize;
} GtMaxbucketinfo;

struct GtLeftborder
{
  uint32_t *uintbounds;
  unsigned long *ulongbounds;
};

struct GtBcktab
{
  GtLeftborder leftborder;
  uint32_t *uintcountspecialcodes;
  unsigned long *ulongcountspecialcodes,
                **distpfxidx,
                sizeofrep;
  unsigned int prefixlength,
               optimalnumofbits;
  GtCodetype numofallcodes,
             numofspecialcodes,
             **multimappower,
             *basepower,
             *filltable;
  GtUchar *qgrambuffer;
  GtMaxbucketinfo maxbucketinfo;
  bool allocated,
       withspecialsuffixes,
       useulong;
  void *mappedptr;
};

void gt_bcktab_leftborder_addcode(GtLeftborder *lb,GtCodetype code)
{
  gt_assert(lb != NULL);
  if (lb->ulongbounds != NULL)
  {
    lb->ulongbounds[code]++;
  } else
  {
    lb->uintbounds[code]++;
  }
}

unsigned long gt_bcktab_leftborder_insertionindex(GtLeftborder *lb,
                                                  GtCodetype code)
{
  gt_assert(lb != NULL);
  if (lb->ulongbounds != NULL)
  {
    return --lb->ulongbounds[code];
  }
  return (unsigned long) --lb->uintbounds[code];
}

void gt_bcktab_leftborder_assign(GtLeftborder *lb,GtCodetype code,
                                 unsigned long value)
{
  gt_assert(lb != NULL);
  if (lb->ulongbounds != NULL)
  {
    lb->ulongbounds[code] = value;
  } else
  {
    gt_assert(value <= (unsigned long) UINT_MAX);
    lb->uintbounds[code] = (uint32_t) value;
  }
}

unsigned long gt_bcktab_get(const GtBcktab *bcktab,GtCodetype code)
{
  gt_assert(bcktab != NULL);
  if (bcktab->leftborder.ulongbounds != NULL)
  {
    return bcktab->leftborder.ulongbounds[code];
  }
  return (unsigned long) bcktab->leftborder.uintbounds[code];
}

static unsigned long numofdistpfxidxcounters(const GtCodetype *basepower,
                                             unsigned int prefixlength)
{
  if (prefixlength > 2U)
  {
    unsigned long numofcounters = 0;
    unsigned int idx;

    for (idx=1U; idx < prefixlength-1; idx++)
    {
      numofcounters += basepower[idx];
    }
    return numofcounters;
  }
  return 0;
}

static uint64_t gt_bcktab_sizeoftable_generic(unsigned int prefixlength,
                                              uint64_t numofallcodes,
                                              uint64_t numofspecialcodes,
                                              unsigned long maxvalue,
                                              const GtCodetype *basepower,
                                              bool withspecialsuffixes)
{
  uint64_t sizeoftable;
  size_t sizeofbasetype;

  sizeofbasetype = maxvalue <= (unsigned long) UINT_MAX
                                 ? sizeof (uint32_t)
                                 : sizeof (unsigned long);
  sizeoftable = (uint64_t) sizeofbasetype * (numofallcodes + 1);
  if (withspecialsuffixes)
  {
    sizeoftable += sizeofbasetype * numofspecialcodes;
    sizeoftable += (uint64_t) sizeof (unsigned long) *
                              numofdistpfxidxcounters(basepower,prefixlength);
  }
  return sizeoftable;
}

uint64_t gt_bcktab_sizeoftable(unsigned int numofchars,
                               unsigned int prefixlength,
                               unsigned long maxvalue,
                               bool withspecialsuffixes)
{
  uint64_t numofallcodes, numofspecialcodes, sizeofbuckettable;
  GtCodetype *basepower;

  numofallcodes = (uint64_t) pow((double) numofchars,(double) prefixlength);
  if (withspecialsuffixes)
  {

    if (prefixlength >= 2U)
    {
      basepower = gt_initbasepower(numofchars,prefixlength-2);
    } else
    {
      basepower = NULL;
    }
    numofspecialcodes = (uint64_t) pow((double) numofchars,
                                       (double) (prefixlength-1));
  } else
  {
    basepower = NULL;
    numofspecialcodes = 0;
  }
  sizeofbuckettable
    = gt_bcktab_sizeoftable_generic(prefixlength,
                                    numofallcodes,
                                    numofspecialcodes,
                                    maxvalue,
                                    basepower,
                                    withspecialsuffixes);
  gt_free(basepower);
  return sizeofbuckettable;
}

unsigned long gt_bcktab_sizeofworkspace(unsigned int prefixlength)
{
  size_t size = sizeof (GtCodetype) * (prefixlength+1);
  size += sizeof (GtCodetype) * prefixlength;
  size += sizeof (GtBcktab);
  size += prefixlength;
  return (unsigned long) size;
}

static void setdistpfxidxptrs(unsigned long **distpfxidx,
                              unsigned long *ptr,
                              const GtCodetype *basepower,
                              unsigned int prefixlength)
{
  unsigned int idx;

  distpfxidx[0] = ptr;
  for (idx=1U; idx<prefixlength-1; idx++)
  {
    distpfxidx[idx] = distpfxidx[idx-1] + basepower[idx];
  }
}

static unsigned long **allocdistprefixindexcounts(const GtCodetype *basepower,
                                                  unsigned int prefixlength)
{
  if (prefixlength > 2U)
  {
    unsigned long numofcounters;

    numofcounters = numofdistpfxidxcounters(basepower,prefixlength);
    if (numofcounters > 0)
    {
      unsigned long *counters, **distpfxidx;

      distpfxidx = gt_malloc(sizeof (unsigned long *) * (prefixlength-1));
      gt_log_log("sizeof (distpfxidx)=%lu bytes",
                  (unsigned long) sizeof (unsigned long *) * (prefixlength-1));
      counters = gt_malloc(sizeof (*counters) * numofcounters);
      gt_log_log("sizeof (counter)=%lu bytes",
                  (unsigned long) sizeof (*counters) * numofcounters);
      memset(counters,0,(size_t) sizeof (*counters) * numofcounters);
      setdistpfxidxptrs(distpfxidx,counters,basepower,prefixlength);
      return distpfxidx;
    }
  }
  return NULL;
}

static GtBcktab *gt_bcktab_new_withinit(unsigned int numofchars,
                                        unsigned int prefixlength,
                                        unsigned long maxvalue,
                                        bool withspecialsuffixes)
{
  GtBcktab *bcktab;
  uint64_t sizeofrep_uint64_t;

  bcktab = gt_malloc(sizeof *bcktab);
  bcktab->leftborder.ulongbounds = NULL;
  bcktab->leftborder.uintbounds = NULL;
  bcktab->ulongcountspecialcodes = NULL;
  bcktab->uintcountspecialcodes = NULL;
  bcktab->distpfxidx = NULL;
  bcktab->mappedptr = NULL;
  bcktab->prefixlength = prefixlength;
  bcktab->withspecialsuffixes = withspecialsuffixes;
  bcktab->basepower = gt_initbasepower(numofchars,prefixlength);
  bcktab->filltable = gt_initfilltable(numofchars,prefixlength);
  bcktab->numofallcodes = bcktab->basepower[prefixlength];
  bcktab->numofspecialcodes = bcktab->basepower[prefixlength-1];
  bcktab->multimappower = gt_initmultimappower(numofchars,prefixlength);
  bcktab->allocated = false;
  bcktab->useulong = false;
  bcktab->qgrambuffer = gt_malloc(sizeof (*bcktab->qgrambuffer) * prefixlength);
  sizeofrep_uint64_t = gt_bcktab_sizeoftable_generic(
                                          prefixlength,
                                          (uint64_t) bcktab->numofallcodes,
                                          (uint64_t) bcktab->numofspecialcodes,
                                          maxvalue,
                                          bcktab->basepower,
                                          withspecialsuffixes);
  bcktab->sizeofrep = CALLCASTFUNC(uint64_t,unsigned_long,sizeofrep_uint64_t);
  return bcktab;
}

GtBcktab *gt_bcktab_new(unsigned int numofchars,
                        unsigned int prefixlength,
                        unsigned long maxvalue,
                        bool storespecialcodes,
                        bool withspecialsuffixes,
                        GtError *err)
{
  GtBcktab *bcktab;
  bool haserr = false;

  bcktab = gt_bcktab_new_withinit(numofchars,prefixlength,maxvalue,
                                  withspecialsuffixes);
  bcktab->allocated = true;
  if (storespecialcodes && bcktab->numofallcodes > 0 &&
      bcktab->numofallcodes-1 > (unsigned long) MAXCODEVALUE)
  {
    gt_error_set(err,"alphasize^prefixlength-1 = " FormatGtCodetype
                  " does not fit into %u"
                  " bits: choose smaller value for prefixlength",
                  bcktab->numofallcodes-1,
                  CODEBITS);
    haserr = true;
  }
  if (!haserr)
  {
    size_t allocsize_bounds, allocsize_countspecialcodes;
    if (maxvalue <= (unsigned long) UINT_MAX)
    {
      allocsize_bounds = sizeof (*bcktab->leftborder.uintbounds) *
                         (bcktab->numofallcodes+1);
      bcktab->leftborder.uintbounds = gt_malloc(allocsize_bounds);
      memset(bcktab->leftborder.uintbounds,0,allocsize_bounds);
      allocsize_countspecialcodes = sizeof (*bcktab->uintcountspecialcodes) *
                                    bcktab->numofspecialcodes;
      bcktab->uintcountspecialcodes = gt_malloc(allocsize_countspecialcodes);
      memset(bcktab->uintcountspecialcodes,0,allocsize_countspecialcodes);
      bcktab->useulong = false;
    } else
    {
      allocsize_bounds = sizeof (*bcktab->leftborder.ulongbounds) *
                         (bcktab->numofallcodes+1);
      bcktab->leftborder.ulongbounds = gt_malloc(allocsize_bounds);
      memset(bcktab->leftborder.ulongbounds,0,allocsize_bounds);
      allocsize_countspecialcodes = sizeof (*bcktab->ulongcountspecialcodes) *
                                    bcktab->numofspecialcodes;
      bcktab->ulongcountspecialcodes = gt_malloc(allocsize_countspecialcodes);
      memset(bcktab->ulongcountspecialcodes,0,allocsize_countspecialcodes);
      bcktab->useulong = true;
    }
    gt_log_log("sizeof (leftborder)=%lu bytes",
                (unsigned long) allocsize_bounds);
    gt_log_log("sizeof (countspecialcodes)=%lu bytes",
                (unsigned long) allocsize_countspecialcodes);
    bcktab->distpfxidx = allocdistprefixindexcounts(bcktab->basepower,
                                                    prefixlength);
    gt_log_log("sizeof (bcktab)=" Formatuint64_t " bytes",
              PRINTuint64_tcast(gt_bcktab_sizeoftable(numofchars,
                                                      prefixlength,
                                                      maxvalue,
                                                      withspecialsuffixes)));
  }
  if (haserr)
  {
    gt_bcktab_delete(bcktab);
    return NULL;
  }
  return bcktab;
}

static void assignbcktabmapspecification(
                                        GtArrayGtMapspecification *mapspectable,
                                        void *voidinfo,
                                        bool writemode)
{
  GtBcktab *bcktab = (GtBcktab *) voidinfo;
  GtMapspecification *mapspecptr;
  unsigned long numofcounters;

  if (bcktab->useulong)
  {
    NEWMAPSPEC(bcktab->leftborder.ulongbounds,GtUlong,
               (unsigned long) (bcktab->numofallcodes+1));
  } else
  {
    NEWMAPSPEC(bcktab->leftborder.uintbounds,Uint32,
               (unsigned long) (bcktab->numofallcodes+1));
  }
  if (bcktab->withspecialsuffixes)
  {
    if (bcktab->useulong)
    {
      NEWMAPSPEC(bcktab->ulongcountspecialcodes,GtUlong,
                 (unsigned long) bcktab->numofspecialcodes);
    } else
    {
      NEWMAPSPEC(bcktab->uintcountspecialcodes,Uint32,
                 (unsigned long) bcktab->numofspecialcodes);
    }
    numofcounters
      = numofdistpfxidxcounters((const GtCodetype *) bcktab->basepower,
                                bcktab->prefixlength);
    if (numofcounters > 0)
    {
      if (!writemode)
      {
        bcktab->distpfxidx = gt_malloc(sizeof (*bcktab->distpfxidx) *
                                       (bcktab->prefixlength-1));
      }
      NEWMAPSPEC(bcktab->distpfxidx[0],GtUlong,numofcounters);
    }
  }
}

int gt_bcktab_flush_to_file(FILE *fp,const GtBcktab *bcktab,GtError *err)
{
  gt_error_check(err);
  return gt_mapspec_flushtheindex2file(fp,
                            assignbcktabmapspecification,
                            (GtBcktab *) bcktab,
                            bcktab->sizeofrep,
                            err);
}

static int fillbcktabmapspecstartptr(GtBcktab *bcktab,
                                     const char *indexname,
                                     GtError *err)
{
  bool haserr = false;
  GtStr *tmpfilename;

  gt_error_check(err);
  tmpfilename = gt_str_new_cstr(indexname);
  gt_str_append_cstr(tmpfilename,BCKTABSUFFIX);
  if (gt_mapspec_fillmapspecstartptr(assignbcktabmapspecification,
                          &bcktab->mappedptr,
                          bcktab,
                          tmpfilename,
                          bcktab->sizeofrep,
                          err) != 0)
  {
    haserr = true;
  }
  gt_str_delete(tmpfilename);
  return haserr ? -1 : 0;
}

GtBcktab *gt_bcktab_map(const char *indexname,
                        unsigned int numofchars,
                        unsigned int prefixlength,
                        unsigned long maxvalue,
                        bool withspecialsuffixes,
                        GtError *err)
{
  GtBcktab *bcktab;

  bcktab = gt_bcktab_new_withinit(numofchars,prefixlength,maxvalue,
                                  withspecialsuffixes);
  bcktab->allocated = false;
  if (fillbcktabmapspecstartptr(bcktab,
                                indexname,
                                err) != 0)
  {
    gt_bcktab_delete(bcktab);
    return NULL;
  }
  if (bcktab->distpfxidx != NULL)
  {
    setdistpfxidxptrs(bcktab->distpfxidx,bcktab->distpfxidx[0],
                      bcktab->basepower,bcktab->prefixlength);
  }
#ifdef SKDEBUG
  gt_bcktab_checkcountspecialcodes(bcktab);
#endif
  return bcktab;
}

void gt_showbcktab(const GtBcktab *bcktab)
{
  unsigned int prefixindex;

  for (prefixindex=1U; prefixindex < bcktab->prefixlength-1; prefixindex++)
  {
    GtCodetype code;
    unsigned long sum = 0;
    for (code = 0; code < bcktab->basepower[prefixindex]; code++)
    {
      sum += bcktab->distpfxidx[prefixindex-1][code];
      printf("distpfxidx[%u][%lu]=%lu\n",
              prefixindex,code,bcktab->distpfxidx[prefixindex-1][code]);
    }
    printf("sum %lu\n",sum);
  }
}

void gt_bcktab_showleftborder(const GtBcktab *bcktab)
{
  GtCodetype idx;

  for (idx=0; idx<bcktab->numofallcodes; idx++)
  {
    printf("leftborder[" FormatGtCodetype "]=%lu\n",idx,
                                                    gt_bcktab_get(bcktab,idx));
  }
}

void gt_bcktab_delete(GtBcktab *bcktab)
{
  if (bcktab == NULL)
  {
    return;
  }
  /*showbcktab(bcktabptr);*/
  if (bcktab->allocated)
  {
    gt_free(bcktab->leftborder.ulongbounds);
    bcktab->leftborder.ulongbounds = NULL;
    gt_free(bcktab->leftborder.uintbounds);
    bcktab->leftborder.uintbounds = NULL;
    gt_free(bcktab->ulongcountspecialcodes);
    gt_free(bcktab->uintcountspecialcodes);
    bcktab->ulongcountspecialcodes = NULL;
    bcktab->uintcountspecialcodes = NULL;
    if (bcktab->distpfxidx != NULL)
    {
      gt_free(bcktab->distpfxidx[0]);
    }
  } else
  {
    if (bcktab->mappedptr != NULL)
    {
      gt_fa_xmunmap(bcktab->mappedptr);
    }
    bcktab->mappedptr = NULL;
    bcktab->leftborder.ulongbounds = NULL;
    bcktab->leftborder.uintbounds = NULL;
    bcktab->ulongcountspecialcodes = NULL;
    bcktab->uintcountspecialcodes = NULL;
    if (bcktab->distpfxidx != NULL)
    {
      bcktab->distpfxidx[0] = NULL;
    }
  }
  gt_free(bcktab->distpfxidx);
  bcktab->distpfxidx = NULL;
  if (bcktab->multimappower != NULL)
  {
    gt_multimappowerfree(&bcktab->multimappower);
  }
  gt_free(bcktab->filltable);
  bcktab->filltable = NULL;
  gt_free(bcktab->basepower);
  bcktab->basepower = NULL;
  gt_free(bcktab->qgrambuffer);
  bcktab->qgrambuffer = NULL;
  gt_free(bcktab);
}

void gt_bcktab_updatespecials(GtBcktab *bcktab,
                              GtCodetype code,
                              unsigned int numofchars,
                              unsigned int prefixindex)
{
  gt_assert(prefixindex > 0);
  gt_assert(prefixindex <= bcktab->prefixlength);
  gt_assert(code < bcktab->numofallcodes);
  if (prefixindex < bcktab->prefixlength-1)
  {
    GtCodetype ordercode = (code - bcktab->filltable[prefixindex])/
                           (bcktab->filltable[prefixindex]+1);
    bcktab->distpfxidx[prefixindex-1][ordercode]++;
  }
  if (bcktab->useulong)
  {
    bcktab->ulongcountspecialcodes[FROMCODE2SPECIALCODE(code,numofchars)]++;
  } else
  {
    bcktab->uintcountspecialcodes[FROMCODE2SPECIALCODE(code,numofchars)]++;
  }
}

void gt_bcktab_addfinalspecials(GtBcktab *bcktab,unsigned int numofchars,
                                unsigned long specialcharacters)
{
  GtCodetype specialcode;

  specialcode = FROMCODE2SPECIALCODE(bcktab->filltable[0],numofchars);
  if (bcktab->useulong)
  {
    bcktab->ulongcountspecialcodes[specialcode]
      += (unsigned long) (specialcharacters + 1);
  } else
  {
    bcktab->uintcountspecialcodes[specialcode]
      += (uint32_t) (specialcharacters + 1);
  }
}

#ifdef SKDEBUG
static unsigned long fromcode2countspecialcodes(GtCodetype code,
                                                const GtBcktab *bcktab)
{
  if (code >= bcktab->filltable[bcktab->prefixlength-1])
  {
    GtCodetype ordercode = code - bcktab->filltable[bcktab->prefixlength-1];
    GtCodetype divisor = bcktab->filltable[bcktab->prefixlength-1] + 1;
    if (ordercode % divisor == 0)
    {
      ordercode /= divisor;
      return bcktab->useulong
               ? bcktab->ulongcountspecialcodes[ordercode]
               : bcktab->uintcountspecialcodes[ordercode];
    }
  }
  return 0;
}

static void pfxidxpartialsums(unsigned long *count,
                              GtCodetype code,
                              const GtBcktab *bcktab)
{
  unsigned int prefixindex;
  unsigned long sum = 0, specialsinbucket;
  GtCodetype ordercode, divisor;

  memset(count,0,sizeof (*count) * (size_t) bcktab->prefixlength);
  for (prefixindex=bcktab->prefixlength-2; prefixindex>=1U; prefixindex--)
  {
    if (code >= bcktab->filltable[prefixindex])
    {
      ordercode = code - bcktab->filltable[prefixindex];
      divisor = bcktab->filltable[prefixindex] + 1;
      if (ordercode % divisor == 0)
      {
        ordercode /= divisor;
        count[prefixindex] = bcktab->distpfxidx[prefixindex-1][ordercode];
        sum += count[prefixindex];
      }
    } else
    {
      break;
    }
  }
  specialsinbucket = fromcode2countspecialcodes(code,bcktab);
  gt_assert(sum <= specialsinbucket);
  count[bcktab->prefixlength-1] = specialsinbucket - sum;
  if (bcktab->prefixlength > 2U)
  {
    for (prefixindex = bcktab->prefixlength-2; prefixindex>=1U; prefixindex--)
    {
      count[prefixindex] += count[prefixindex+1];
    }
#ifndef NDEBUG
    if (specialsinbucket != count[1])
    {
      fprintf(stderr,"code " FormatGtCodetype ": sum = %lu != %lu = count[1]\n",
              code,sum,count[1]);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
#endif
  }
}

void checkcountspecialcodes(const GtBcktab *bcktab)
{
  GtCodetype code;
  unsigned long *count;

  if (bcktab->prefixlength >= 2U)
  {
    count = gt_malloc(sizeof (*count) * bcktab->prefixlength);
    for (code=0; code<bcktab->numofallcodes; code++)
    {
      pfxidxpartialsums(count,
                        code,
                        bcktab);
    }
    gt_free(count);
  }
}
#endif

GtCodetype gt_bcktab_codedownscale(const GtBcktab *bcktab,
                                   GtCodetype code,
                                   unsigned int prefixindex,
                                   unsigned int maxprefixlen)
{
  unsigned int remain;

  code -= bcktab->filltable[maxprefixlen];
  remain = maxprefixlen-prefixindex;
  code %= (bcktab->filltable[remain]+1);
  code *= bcktab->basepower[remain];
  code += bcktab->filltable[prefixindex];
  return code;
}

unsigned int gt_bcktab_calcboundsparts(GtBucketspecification *bucketspec,
                                       const GtBcktab *bcktab,
                                       GtCodetype code,
                                       GtCodetype maxcode,
                                       unsigned long totalwidth,
                                       unsigned int rightchar,
                                       unsigned int numofchars)
{
  bucketspec->left = gt_bcktab_get(bcktab,code);
  if (code == maxcode)
  {
    gt_assert(totalwidth >= bucketspec->left);
    bucketspec->nonspecialsinbucket
      = (unsigned long) (totalwidth - bucketspec->left);
  } else
  {
    if (gt_bcktab_get(bcktab,code+1) > 0)
    {
      bucketspec->nonspecialsinbucket
        = gt_bcktab_get(bcktab,code+1) - bucketspec->left;
    } else
    {
      bucketspec->nonspecialsinbucket = 0;
    }
  }
  gt_assert(rightchar == (unsigned int) (code % numofchars));
  if (rightchar == numofchars - 1)
  {
    bucketspec->specialsinbucket
      = bcktab->useulong
         ? bcktab->ulongcountspecialcodes[FROMCODE2SPECIALCODE(code,numofchars)]
         : (unsigned long) bcktab->uintcountspecialcodes[
                                   FROMCODE2SPECIALCODE(code,numofchars)];
    if (bucketspec->nonspecialsinbucket >= bucketspec->specialsinbucket)
    {
      bucketspec->nonspecialsinbucket -= bucketspec->specialsinbucket;
    } else
    {
      bucketspec->nonspecialsinbucket = 0;
    }
    rightchar = 0;
  } else
  {
    bucketspec->specialsinbucket = 0;
    rightchar++;
  }
  return rightchar;
}

void gt_bcktab_calcboundaries(GtBucketspecification *bucketspec,
                              const GtBcktab *bcktab,
                              GtCodetype code)
{
  unsigned int numofchars = (unsigned int) bcktab->basepower[1];
  gt_assert(code != bcktab->numofallcodes);
  (void) gt_bcktab_calcboundsparts(bucketspec,
                                   bcktab,
                                   code,
                                   bcktab->numofallcodes, /*code<numofallcodes*/
                                   0,  /* not necessary */
                                   (unsigned int) (code % numofchars),
                                   numofchars);
}

unsigned long gt_bcktab_calcrightbounds(const GtBcktab *bcktab,
                                        GtCodetype code,
                                        GtCodetype maxcode,
                                        unsigned long totalwidth)
{
  return (code == maxcode) ? totalwidth : gt_bcktab_get(bcktab,code+1);
}

void gt_bcktab_determinemaxsize(GtBcktab *bcktab,
                                const GtCodetype mincode,
                                const GtCodetype maxcode,
                                unsigned long partwidth,
                                unsigned int numofchars)
{
  unsigned int rightchar = (unsigned int) (mincode % numofchars);
  GtBucketspecification bucketspec;
  GtCodetype code;

#ifdef SKDEBUG
  printf("mincode=%lu,maxcode=%lu,partwidth=%lu,totallength=%lu\n",
          (unsigned long) mincode,(unsigned long) maxcode,
          partwidth,totallength);
#endif
  bcktab->maxbucketinfo.specialsmaxbucketsize = 1UL;
  bcktab->maxbucketinfo.nonspecialsmaxbucketsize = 1UL;
  bcktab->maxbucketinfo.maxbucketsize = 1UL;
  for (code = mincode; code <= maxcode; code++)
  {
    rightchar = gt_bcktab_calcboundsparts(&bucketspec,
                                          bcktab,
                                          code,
                                          maxcode,
                                          partwidth,
                                          rightchar,
                                          numofchars);
    if (bucketspec.nonspecialsinbucket >
        bcktab->maxbucketinfo.nonspecialsmaxbucketsize)
    {
      bcktab->maxbucketinfo.nonspecialsmaxbucketsize
        = bucketspec.nonspecialsinbucket;
    }
    if (bucketspec.specialsinbucket >
        bcktab->maxbucketinfo.specialsmaxbucketsize)
    {
      bcktab->maxbucketinfo.specialsmaxbucketsize = bucketspec.specialsinbucket;
    }
    if (bucketspec.nonspecialsinbucket + bucketspec.specialsinbucket
        > bcktab->maxbucketinfo.maxbucketsize)
    {
      bcktab->maxbucketinfo.maxbucketsize = bucketspec.nonspecialsinbucket +
                                            bucketspec.specialsinbucket;
    }
  }
  /*
  gt_logger_log(logger,"maxbucket (specials)=%lu",
              bcktab->maxbucketinfo.specialsmaxbucketsize);
  gt_logger_log(logger,"maxbucket (nonspecials)=%lu",
              bcktab->maxbucketinfo.nonspecialsmaxbucketsize);
  gt_logger_log(logger,"maxbucket (all)=%lu",
              bcktab->maxbucketinfo.maxbucketsize);
  */
}

static unsigned long gt_bcktab_nonspecialsmaxsize(const GtBcktab *bcktab)
{
  return bcktab->maxbucketinfo.nonspecialsmaxbucketsize;
}

unsigned int gt_bcktab_prefixlength(const GtBcktab *bcktab)
{
  return bcktab->prefixlength;
}

unsigned int gt_bcktab_singletonmaxprefixindex(const GtBcktab *bcktab,
                                               GtCodetype code)
{
  if (bcktab->prefixlength > 2U)
  {
    GtCodetype ordercode, divisor;
    unsigned int prefixindex;

    for (prefixindex=bcktab->prefixlength-2; prefixindex>=1U; prefixindex--)
    {
      if (code >= bcktab->filltable[prefixindex])
      {
        ordercode = code - bcktab->filltable[prefixindex];
        divisor = bcktab->filltable[prefixindex] + 1;
        if (ordercode % divisor == 0)
        {
          ordercode /= divisor;
          if (bcktab->distpfxidx[prefixindex-1][ordercode] > 0)
          {
            return prefixindex;
          }
        }
      } else
      {
        break;
      }
    }
  }
  return bcktab->prefixlength-1;
}

/* The following function is not used */

unsigned long gt_bcktab_distpfxidxpartialsums(const GtBcktab *bcktab,
                                              GtCodetype code,
                                              unsigned int lowerbound)
{
  GtCodetype ordercode, divisor;
  unsigned long sum = 0;
  unsigned int prefixindex;

  for (prefixindex=bcktab->prefixlength-2; prefixindex>lowerbound;
       prefixindex--)
  {
    if (code >= bcktab->filltable[prefixindex])
    {
      ordercode = code - bcktab->filltable[prefixindex];
      divisor = bcktab->filltable[prefixindex] + 1;
      if (ordercode % divisor == 0)
      {
        ordercode /= divisor;
        sum += bcktab->distpfxidx[prefixindex-1][ordercode];
      }
    } else
    {
      break;
    }
  }
  return sum;
}

unsigned int gt_bcktab_pfxidx2lcpvalues_uint8(unsigned int *minprefixindex,
                                              uint8_t *smalllcpvalues,
                                              unsigned long specialsinbucket,
                                              const GtBcktab *bcktab,
                                              GtCodetype code)
{
  unsigned int prefixindex, maxprefixindex = 0;
  unsigned long idx, insertpos;
  GtCodetype ordercode, divisor;

  gt_assert(smalllcpvalues != NULL);
  *minprefixindex = bcktab->prefixlength;
  insertpos = specialsinbucket;
  for (prefixindex=1U; prefixindex<bcktab->prefixlength-1; prefixindex++)
  {
    if (code >= bcktab->filltable[prefixindex])
    {
      ordercode = code - bcktab->filltable[prefixindex];
      divisor = bcktab->filltable[prefixindex] + 1;
      if (ordercode % divisor == 0)
      {
        ordercode /= divisor;
        if (bcktab->distpfxidx[prefixindex-1][ordercode] > 0)
        {
          maxprefixindex = prefixindex;
          if (*minprefixindex > prefixindex)
          {
            *minprefixindex = prefixindex;
          }
          for (idx=0; idx < bcktab->distpfxidx[prefixindex-1][ordercode]; idx++)
          {
            gt_assert(insertpos > 0);
            smalllcpvalues[--insertpos] = (uint8_t) prefixindex;
          }
        }
      }
    }
  }
  if (insertpos > 0)
  {
    maxprefixindex = bcktab->prefixlength-1;
    if (*minprefixindex == bcktab->prefixlength)
    {
      *minprefixindex = bcktab->prefixlength-1;
    }
    while (insertpos > 0)
    {
      smalllcpvalues[--insertpos] = (uint8_t) (bcktab->prefixlength-1);
    }
  }
  return maxprefixindex;
}

unsigned int gt_bcktab_pfxidx2lcpvalues_ulong(unsigned int *minprefixindex,
                                              unsigned long *bucketoflcpvalues,
                                              unsigned long specialsinbucket,
                                              const GtBcktab *bcktab,
                                              GtCodetype code)
{
  unsigned int prefixindex, maxprefixindex = 0;
  unsigned long idx, insertpos;
  GtCodetype ordercode, divisor;

  gt_assert(bucketoflcpvalues != NULL);
  *minprefixindex = bcktab->prefixlength;
  insertpos = specialsinbucket;
  for (prefixindex=1U; prefixindex<bcktab->prefixlength-1; prefixindex++)
  {
    if (code >= bcktab->filltable[prefixindex])
    {
      ordercode = code - bcktab->filltable[prefixindex];
      divisor = bcktab->filltable[prefixindex] + 1;
      if (ordercode % divisor == 0)
      {
        ordercode /= divisor;
        if (bcktab->distpfxidx[prefixindex-1][ordercode] > 0)
        {
          maxprefixindex = prefixindex;
          if (*minprefixindex > prefixindex)
          {
            *minprefixindex = prefixindex;
          }
          for (idx=0; idx < bcktab->distpfxidx[prefixindex-1][ordercode]; idx++)
          {
            gt_assert(insertpos > 0);
            bucketoflcpvalues[--insertpos] = (unsigned long) prefixindex;
          }
        }
      }
    }
  }
  if (insertpos > 0)
  {
    maxprefixindex = bcktab->prefixlength-1;
    if (*minprefixindex == bcktab->prefixlength)
    {
      *minprefixindex = bcktab->prefixlength-1;
    }
    while (insertpos > 0)
    {
      bucketoflcpvalues[--insertpos] = (unsigned long) (bcktab->prefixlength-1);
    }
  }
  return maxprefixindex;
}

const GtCodetype **gt_bcktab_multimappower(const GtBcktab *bcktab)
{
  return (const GtCodetype **) bcktab->multimappower;
}

GtCodetype gt_bcktab_filltable(const GtBcktab *bcktab,unsigned int idx)
{
  return bcktab->filltable[idx];
}

GtLeftborder *gt_bcktab_leftborder(GtBcktab *bcktab)
{
  return &bcktab->leftborder;
}

GtCodetype gt_bcktab_numofallcodes(const GtBcktab *bcktab)
{
  return bcktab->numofallcodes;
}

unsigned long gt_bcktab_leftborderpartialsums(
                             unsigned long *saved_bucketswithoutwholeleaf,
                             unsigned long *numofsuffixestosort,
                             GtBcktab *bcktab,
                             const GtBitsequence *markwholeleafbuckets)
{
  unsigned long code, largestbucketsize, sumbuckets, saved = 0, currentsize;

  gt_assert(bcktab->numofallcodes > 0);
  if (markwholeleafbuckets != NULL && !GT_ISIBITSET(markwholeleafbuckets,0))
  {
    saved += gt_bcktab_get(bcktab,0);
    gt_bcktab_leftborder_assign(&bcktab->leftborder,0,0);
  }
  largestbucketsize = sumbuckets = gt_bcktab_get(bcktab,0);
  for (code = 1UL; code < bcktab->numofallcodes; code++)
  {
    currentsize = gt_bcktab_get(bcktab,code);
    if (markwholeleafbuckets != NULL &&
        !GT_ISIBITSET(markwholeleafbuckets,code))
    {
      saved += currentsize;
      gt_bcktab_leftborder_assign(&bcktab->leftborder,code,
                                  gt_bcktab_get(bcktab,code-1));
    } else
    {
      sumbuckets += currentsize;
      if (largestbucketsize < currentsize)
      {
        largestbucketsize = currentsize;
      }
      gt_bcktab_leftborder_assign(&bcktab->leftborder,code,
                                  gt_bcktab_get(bcktab,code) +
                                  gt_bcktab_get(bcktab,code-1));
    }
  }
  gt_bcktab_leftborder_assign(&bcktab->leftborder,bcktab->numofallcodes,
                              sumbuckets);
  if (saved_bucketswithoutwholeleaf != NULL)
  {
    *saved_bucketswithoutwholeleaf = saved;
  }
  if (numofsuffixestosort != NULL)
  {
    *numofsuffixestosort = sumbuckets;
  }
  return largestbucketsize;
}

size_t gt_bcktab_sizeforlcpvalues(const GtBcktab *bcktab)
{
  size_t sizespeciallcps, sizelcps;

  sizespeciallcps = sizeof (uint8_t) *
                    bcktab->maxbucketinfo.specialsmaxbucketsize;
  sizelcps
    = sizeof (unsigned long) * gt_bcktab_nonspecialsmaxsize(bcktab);
  return MAX(sizelcps,sizespeciallcps);
}

GtCodetype gt_bcktab_findfirstlarger(const GtBcktab *bcktab,
                                     unsigned long suftaboffset)
{
  GtCodetype left = 0, right = bcktab->numofallcodes, mid,
             found = bcktab->numofallcodes;
  unsigned long midval;

  while (left+1 < right)
  {
    mid = GT_DIV2(left+right);
    midval = gt_bcktab_get(bcktab,mid);
    if (suftaboffset == midval)
    {
      return mid;
    }
    if (suftaboffset < midval)
    {
      found = mid;
      right = mid - 1;
    } else
    {
      left = mid + 1;
    }
  }
  return found;
}

#ifdef SKDEBUG
#include "qgram2code.h"
void gt_bcktab_consistencyofsuffix(int line,
                                   const GtEncseq *encseq,
                                   GtReadmode readmode,
                                   const GtBcktab *bcktab,
                                   unsigned int numofchars,
                                   const Suffixwithcode *suffix)
{
  unsigned int idx, firstspecial = bcktab->prefixlength, gramfirstspecial;
  GtCodetype qgramcode = 0;
  unsigned long totallength;
  GtUchar cc = 0;

  totallength = gt_encseq_total_length(encseq);
  for (idx=0; idx<bcktab->prefixlength; idx++)
  {
    if (suffix->startpos + idx >= totallength)
    {
      firstspecial = idx;
      break;
    }
    cc = gt_encseq_get_encoded_char(encseq,suffix->startpos + idx,
                                           readmode);
    if (ISSPECIAL(cc))
    {
      firstspecial = idx;
      break;
    }
    bcktab->qgrambuffer[idx] = cc;
  }
  for (idx=firstspecial; idx<bcktab->prefixlength; idx++)
  {
    bcktab->qgrambuffer[idx] = (GtUchar) (numofchars-1);
  }
  gramfirstspecial = qgram2code(&qgramcode,
                                (const GtCodetype **) bcktab->multimappower,
                                bcktab->prefixlength,
                                bcktab->qgrambuffer);
  gt_assert(gramfirstspecial == bcktab->prefixlength);
  gt_assert(qgramcode == suffix->code);
  if (firstspecial != suffix->prefixindex)
  {
    fprintf(stderr,"line %d: code=%u: ",line,suffix->code);
    fprintf(stderr,"firstspecial = %u != %u = suffix->prefixindex\n",
                    firstspecial,suffix->prefixindex);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}
#endif
