/*
  Copyright (c) 2007      Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>,
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include <cairo.h>
#include <cairo-pdf.h>
#include <cairo-ps.h>
#ifdef CAIRO_HAS_SVG_SURFACE
#include <cairo-svg.h>
#endif
#include "core/fileutils.h"
#include "core/genfile.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/str.h"
#include "annotationsketch/graphics.h"

struct Graphics {
  cairo_t *cr;
  cairo_surface_t *surf;
  Str *outbuf;
  GraphicsOutType type;
  double margin_x, margin_y, height, width;
};

static cairo_status_t str_write_func(void *closure, const unsigned char *data,
                                     unsigned int length)
{
  Str *stream = closure;
  assert(stream);
  str_append_cstr_nt(stream, (char*) data, length);
  return CAIRO_STATUS_SUCCESS;
}

void graphics_initialize(Graphics *g, GraphicsOutType type,
                         unsigned int width, unsigned int height)
{
  g->outbuf = str_new();
  switch (type)
  {
    case GRAPHICS_PDF:
#ifdef CAIRO_HAS_PDF_SURFACE
      g->surf = cairo_pdf_surface_create_for_stream(str_write_func,
                                                    g->outbuf,
                                                    width,
                                                    height);
#endif
      break;
    case GRAPHICS_PS:
#ifdef CAIRO_HAS_PS_SURFACE
      g->surf = cairo_ps_surface_create_for_stream(str_write_func,
                                                   g->outbuf,
                                                   width,
                                                   height);
#endif
      break;
    case GRAPHICS_SVG:
#ifdef CAIRO_HAS_SVG_SURFACE
      g->surf = cairo_svg_surface_create_for_stream(str_write_func,
                                                    g->outbuf,
                                                    width,
                                                    height);
#endif
      break;
    case GRAPHICS_PNG:
    default:
      g->surf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
      break;
  }
  assert(g->surf && cairo_surface_status(g->surf) == CAIRO_STATUS_SUCCESS);
  g->cr = cairo_create(g->surf);
  assert(cairo_status(g->cr) == CAIRO_STATUS_SUCCESS);
  g->width = width;
  g->height = height;
  g->margin_x = g->margin_y = 20;
  cairo_set_source_rgba(g->cr, 1, 1, 1, 1);
  cairo_paint(g->cr);
  cairo_set_line_join(g->cr, CAIRO_LINE_JOIN_ROUND);
  cairo_set_line_cap(g->cr, CAIRO_LINE_CAP_ROUND);
  cairo_select_font_face(g->cr, "sans", CAIRO_FONT_SLANT_NORMAL,
                         CAIRO_FONT_WEIGHT_NORMAL);
  g->type = type;
}

Graphics* graphics_new(GraphicsOutType type,
                       unsigned int width, unsigned int height)
{
  Graphics *g = ma_calloc(1, sizeof (Graphics));
  graphics_initialize(g, type, width, height);
  return g;
}

void graphics_set_font(Graphics *g, const char *family,
                       FontSlant slant, FontWeight weight)
{
  assert(g && family);
  cairo_select_font_face(g->cr, family,
                         (cairo_font_slant_t) slant,
                         (cairo_font_weight_t) weight);
}

void graphics_draw_text(Graphics *g, double x, double y, const char *text)
{
  cairo_text_extents_t ext;
  assert(g && text);
  cairo_text_extents(g->cr, text, &ext);
  if (x+ext.width > g->width)
    return;
  cairo_set_source_rgb(g->cr, 0, 0, 0);
  cairo_move_to(g->cr, x, y);
  cairo_show_text(g->cr, text);
}

void graphics_draw_text_centered(Graphics *g, double x, double y,
                                 const char *text)
{
  cairo_text_extents_t ext;
  assert(g && text);
  cairo_set_source_rgb(g->cr, 0, 0, 0);
  /* get text extents */
  cairo_text_extents(g->cr, text, &ext);
  /* draw text w/ its center at the given coords */
  cairo_move_to(g->cr, x-(ext.width/2)-1, y);
  cairo_show_text(g->cr, text);
}

void graphics_draw_text_right(Graphics *g, double x, double y, const char *text)
{
  cairo_text_extents_t ext;
  assert(g && text);
  cairo_set_source_rgb(g->cr, 0, 0, 0);
  /* get text extents */
  cairo_text_extents(g->cr, text, &ext);
  /* draw text w/ its right end at the given coords */
  cairo_move_to(g->cr, x-(ext.width)-1, y);
  cairo_show_text(g->cr, text);
}

void graphics_draw_colored_text(Graphics *g, double x, double y, GT_Color color,
                                const char *text)
{
  assert(g && text);
  cairo_set_source_rgb(g->cr, color.red, color.green, color.blue);
  cairo_move_to(g->cr, x, y);
  cairo_show_text(g->cr, text);
}

double graphics_get_image_height(Graphics *g)
{
  assert(g);
  return g->height;
}

double graphics_get_image_width(Graphics *g)
{
  assert(g);
  return g->width;
}

void graphics_set_margins(Graphics *g, double margin_x, double margin_y)
{
  assert(g);
  g->margin_x = margin_x;
  g->margin_y = margin_y;
}

void graphics_draw_horizontal_line(Graphics *g, double x, double y,
                                   double width)
{
  assert(g);
  cairo_move_to(g->cr, x, y);
  cairo_rel_line_to(g->cr, width, 0);
  cairo_stroke(g->cr);
}

void graphics_draw_vertical_line(Graphics *g, double x, double y,
                                 GT_Color color, double length)
{
  assert(g);
  cairo_save(g->cr);
  cairo_move_to(g->cr, x, y);
  cairo_set_line_width(g->cr, 0.7);
  cairo_set_source_rgb(g->cr, color.red, color.green, color.blue);
  cairo_rel_line_to(g->cr, 0, length);
  cairo_stroke(g->cr);
  cairo_restore(g->cr);
}

double graphics_get_text_width(Graphics *g, const char* text)
{
  cairo_text_extents_t ext;
  assert(g);
  /* get text extents */
  cairo_text_extents(g->cr, text, &ext);
  return ext.width;
}

double graphics_get_text_height(Graphics *g)
{
  cairo_font_extents_t ext;
  assert(g);
  cairo_font_extents(g->cr, &ext);
  return ext.height;
}

void graphics_draw_box(Graphics *g, double x, double y, double width,
                       double height, GT_Color fill_color,
                       ArrowStatus arrow_status, double arrow_width,
                       double stroke_width, GT_Color stroke_color, bool dashed)
{
  assert(g);
  double dashes[]={2.0};
  /* save cairo context */
  cairo_save(g->cr);
  cairo_rectangle(g->cr, g->margin_x, g->margin_y,
                  g->width-2*g->margin_x, g->height-2*g->margin_y);
  cairo_clip(g->cr);
  /* construct shape of the box or arrow */
  switch (arrow_status)
  {
    case ARROW_RIGHT:
      cairo_move_to(g->cr, x, y);
      if (width - arrow_width > 0)
        cairo_rel_line_to(g->cr, width - arrow_width, 0);
      cairo_line_to(g->cr, x + width, y + height / 2);
      if (width - arrow_width > 0)
        cairo_line_to(g->cr, x + width - arrow_width, y + height);
      cairo_line_to(g->cr, x, y + height);
      cairo_close_path(g->cr);
      break;
    case ARROW_LEFT:
      cairo_move_to(g->cr, x + width, y);
      if (width - arrow_width > 0)
        cairo_rel_line_to(g->cr, -(width - arrow_width), 0);
      cairo_line_to(g->cr, x, y + height / 2);
      cairo_line_to(g->cr, x + MIN(width, arrow_width), y + height);
      if (width - arrow_width > 0)
        cairo_line_to(g->cr, x + width, y + height);
      cairo_close_path(g->cr);
      break;
    case ARROW_BOTH: /* XXX */
    case ARROW_NONE:
      cairo_rectangle(g->cr, x, y, width, height);
   }
   /* fill area */
   cairo_set_source_rgba(g->cr, fill_color.red,
                                fill_color.green,
                                fill_color.blue,
                                0.5);
   cairo_fill_preserve(g->cr);
   /* draw outline */
   cairo_set_line_width(g->cr, stroke_width);
   cairo_set_source_rgba(g->cr, stroke_color.red,
                                stroke_color.green,
                                stroke_color.blue,
                                0.7);
   if (dashed)
     cairo_set_dash(g->cr, dashes, 1, (double) 0);
   cairo_stroke(g->cr);
   /* restore cairo context */
   cairo_restore(g->cr);
}

void graphics_draw_dashes(Graphics *g, double x, double y, double width,
                          double height, ArrowStatus arrow_status,
                          double arrow_width, double stroke_width,
                          GT_Color stroke_color)
{
  double dashes[] = {3.0};
  assert(g);
  /* save cairo context */
  cairo_save(g->cr);
  cairo_rectangle(g->cr, g->margin_x, g->margin_y,
                  g->width-2*g->margin_x, g->height-2*g->margin_y);
  cairo_clip(g->cr);
  cairo_set_line_width(g->cr, stroke_width);
  cairo_set_source_rgb(g->cr, stroke_color.red,
                              stroke_color.green,
                              stroke_color.blue);
  /* create arrowhead path */
  if (width - arrow_width > 0)
  {
    switch (arrow_status)
    {
      case ARROW_RIGHT:
        cairo_move_to(g->cr, x + width - arrow_width, y);
        cairo_line_to(g->cr, x + width, y + height / 2);
        cairo_line_to(g->cr, x + width - arrow_width, y + height);
        /* draw arrowhead */
        cairo_stroke(g->cr);
        break;
      case ARROW_LEFT:
        cairo_move_to(g->cr, width - arrow_width, y);
        cairo_line_to(g->cr, x + width, y + height / 2);
        cairo_line_to(g->cr, x + width - arrow_width, y + height);
        /* draw arrowhead */
        cairo_stroke(g->cr);
        break;
      case ARROW_BOTH: /* XXX */
      case ARROW_NONE: break;
    }
  }
  /* draw dashes */
  cairo_set_dash(g->cr, dashes, 1, (double) 0);
  cairo_move_to(g->cr, x, y+(height/2));
  cairo_rel_line_to(g->cr, width, 0);
  cairo_stroke(g->cr);
  /* restore cairo context */
  cairo_restore(g->cr);
}

void graphics_draw_caret(Graphics *g, double x, double y, double width,
                         double height, ArrowStatus arrow_status,
                         double arrow_width, double stroke_width,
                         GT_Color stroke_color)
{
  assert(g);
  /* save cairo context */
  cairo_save(g->cr);
  cairo_rectangle(g->cr, g->margin_x, g->margin_y,
                  g->width-2*g->margin_x, g->height-2*g->margin_y);
  cairo_clip(g->cr);
  /* create caret path */
  switch (arrow_status)
  {
    case ARROW_RIGHT:
      if (width - arrow_width > 0)
      {
        width -= arrow_width;
        cairo_move_to(g->cr, x, y + (height/2));
        cairo_line_to(g->cr, x + (width/2), y);
        cairo_line_to(g->cr, x + width, y + (height/2));
        cairo_line_to(g->cr, x + width + arrow_width, y + (height/2));
        /* add arrowhead */
        cairo_move_to(g->cr, x + width, y);
        cairo_line_to(g->cr, x + width + arrow_width, y + height / 2);
        cairo_line_to(g->cr, x + width, y + height);
      }
      break;
    case ARROW_LEFT:
      if (width - arrow_width > 0)
      {
        width -= arrow_width;
        /* draw arrowhead */
        cairo_move_to(g->cr, x + arrow_width , y);
        cairo_line_to(g->cr, x, y + (height/2));
        cairo_line_to(g->cr, x + arrow_width, y + height);
        cairo_move_to(g->cr, x, y + (height/2));
        cairo_line_to(g->cr, x + arrow_width, y + (height/2));
        cairo_line_to(g->cr, x + arrow_width + (width/2), y);
        cairo_line_to(g->cr, x + arrow_width + width, y + (height/2));
      }
      break;
    case ARROW_BOTH: /* XXX */
    case ARROW_NONE:
     cairo_move_to(g->cr, x, y + (height/2));
     cairo_line_to(g->cr, x + (width/2), y);
     cairo_line_to(g->cr, x + width, y + (height/2));
   }
   /* draw caret */
   cairo_set_line_width(g->cr, stroke_width);
   cairo_set_source_rgb(g->cr, stroke_color.red,
                               stroke_color.green,
                               stroke_color.blue);
   cairo_stroke(g->cr);
   /* restore cairo context */
   cairo_restore(g->cr);
}

void graphics_draw_rectangle(Graphics *g, double x, double y,
                             bool filled, GT_Color fill_color, bool outlined,
                             GT_Color outline_color, double outline_width,
                             double width)
{
  assert(g);
  /* save cairo context */
  cairo_save(g->cr);
  cairo_new_path(g->cr);
  cairo_rectangle(g->cr, x, y, x + width, y + width);
  if (filled) {
    cairo_set_source_rgb(g->cr, fill_color.red,
                                fill_color.green,
                                fill_color.blue);
    cairo_fill_preserve(g->cr);
  }
  if (outlined) {
    cairo_set_line_width(g->cr, outline_width);
    cairo_set_source_rgb(g->cr, outline_color.red,
                                outline_color.green,
                                outline_color.blue);
    cairo_stroke(g->cr);
  }
  /* restore cairo context */
  cairo_restore(g->cr);
}

void graphics_draw_arrowhead(Graphics *g, double x, double y,
                             GT_Color color, ArrowStatus arrow_status)
{
  double arrow_height = 8, arrow_width = 5;
  assert(g);
  /* save cairo context */
  cairo_save(g->cr);
  cairo_reset_clip(g->cr);
  cairo_set_source_rgb(g->cr, color.red, color.green, color.blue);
  switch (arrow_status)
  {
    case ARROW_LEFT:
      cairo_move_to(g->cr, x+arrow_width, y);
      cairo_line_to(g->cr, x, y+(arrow_height/2));
      cairo_line_to(g->cr, x+arrow_width, y+arrow_height);
      cairo_close_path(g->cr);
      /* fill area */
      cairo_fill_preserve(g->cr);
      cairo_stroke(g->cr);
      break;
    case ARROW_RIGHT:
      cairo_move_to(g->cr, x, y);
      cairo_line_to(g->cr, x+arrow_width, y+(arrow_height/2));
      cairo_line_to(g->cr, x, y+arrow_height);
      cairo_close_path(g->cr);
      /* fill area */
      cairo_fill_preserve(g->cr);
      cairo_stroke(g->cr);
      break;
    case ARROW_BOTH: /* XXX */
    case ARROW_NONE: break;
  }
  /* restore cairo context */
  cairo_restore(g->cr);
}

int graphics_save_to_file(const Graphics *g, const char *filename, Error *err)
{
  cairo_status_t rval;
  error_check(err);
  GenFile *outfile;
  assert(g && filename);

  switch (g->type)
  {
    case GRAPHICS_PNG:
      rval = cairo_surface_write_to_png(g->surf, filename);
      assert(rval == CAIRO_STATUS_SUCCESS || rval == CAIRO_STATUS_WRITE_ERROR);
      if (rval == CAIRO_STATUS_WRITE_ERROR)
      {
        error_set(err, "an I/O error occurred while attempting "
                       "to write image file \"%s\"", filename);
        return -1;
      }
      break;
    default:
      cairo_show_page(g->cr);
      cairo_surface_flush(g->surf);
      cairo_surface_finish(g->surf);
      outfile = genfile_open(GFM_UNCOMPRESSED, filename, "w+", err);
      if (!error_is_set(err))
      {
        genfile_xwrite(outfile, str_get_mem(g->outbuf), str_length(g->outbuf));
        genfile_close(outfile);
      } else return -1;
      break;
  }
  return 0;
}

void graphics_save_to_stream(const Graphics *g, Str *stream)
{
  cairo_status_t rval;
  assert(g && stream);
  rval = cairo_surface_write_to_png_stream(g->surf, str_write_func, stream);
  assert(rval == CAIRO_STATUS_SUCCESS); /* str_write_func() is sane */
}

void graphics_delete(Graphics *g)
{
  if (!g) return;
  cairo_surface_destroy(g->surf); /* reference counted */
  cairo_destroy(g->cr);
  if (g->outbuf)
    str_delete(g->outbuf);
  ma_free(g);
}
