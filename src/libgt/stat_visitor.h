/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef STAT_VISITOR_H
#define STAT_VISITOR_H

/* implements the ``genome visitor'' interface, gathers statistics */
typedef struct StatVisitor StatVisitor;

#include "genome_visitor.h"

const GenomeVisitorClass* stat_visitor_class(void);
GenomeVisitor*            stat_visitor_new(unsigned int gene_length_distri,
                                           unsigned int gene_score_distri);
void                      stat_visitor_show_stats(GenomeVisitor*);

#endif
