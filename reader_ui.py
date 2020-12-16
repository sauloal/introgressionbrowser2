#!/usr/bin/env python3

import os
import sys
import typing
import math
import time

import numpy as np
# from matplotlib.cm import ScalarMappable
from matplotlib.cm import get_cmap, ScalarMappable
from matplotlib.colors import Normalize
#, Colormap

import flexx
from   flexx import flx, ui

import reader

DEBUG       = True


def setHidden(box, status):
    if " hidden" in box.css_class:
        if status:
            pass
        else:
            box.set_css_class(box.css_class.replace(" hidden", ""))

    else:
        if status:
            box.set_css_class(box.css_class + " hidden")
        else:
            pass

color_maps = {
    'Perceptually Uniform Sequential': [
        'viridis', 'plasma', 'inferno', 'magma',
        'cividis'
    ],
    'Sequential': [
        'Greys'  , 'Purples', 'Blues' , 'Greens',
        'Oranges', 'Reds'   , 'YlOrBr', 'YlOrRd',
        'OrRd'   , 'PuRd'   , 'RdPu'  , 'BuPu'  ,
        'GnBu'   , 'PuBu'   , 'YlGnBu', 'PuBuGn',
        'BuGn'   , 'YlGn'
    ],
    'Sequential (2)': [
        'binary', 'gist_yarg', 'gist_gray', 'gray'  ,
        'bone'  , 'pink'     , 'spring'   , 'summer',
        'autumn', 'winter'   , 'cool'     , 'Wistia',
        'hot'   , 'afmhot'   , 'gist_heat', 'copper'
    ],
    'Diverging': [
        'PiYG'    , 'PRGn'    , 'BrBG'  , 'PuOr'  ,
        'RdGy'    , 'RdBu'    , 'RdYlBu', 'RdYlGn',
        'Spectral', 'coolwarm', 'bwr'   , 'seismic'
    ],
    'Cyclic': [
        'twilight', 'twilight_shifted', 'hsv'
    ],
    'Qualitative': [
        'Pastel1', 'Pastel2', 'Paired', 'Accent',
        'Dark2'  , 'Set1'   , 'Set2'  , 'Set3'  ,
        'tab10'  , 'tab20'  , 'tab20b', 'tab20c'
    ],
    'Miscellaneous': [
        'flag'        , 'prism'        , 'ocean'        , 'gist_earth'  ,
        'terrain'     , 'gist_stern'   , 'gnuplot'      , 'gnuplot2'    ,
        'CMRmap'      , 'cubehelix'    , 'brg'          , 'gist_rainbow',
        'rainbow'     , 'jet'          , 'nipy_spectral', 'gist_ncar'
    ]
}

def get_color_maps(map_name, min_val=0.0, max_val=1.0, levels=256, bad='pink', over='pink', under='pink'):
    """
        norm, cmap, f = get_color_maps("Greys", 0, 100)
        print(-1, f(-1.0))
        print(0, f(0.0))
        print(0.5, f(0.5))
        print(1, f(1.0))
        print(10, f(10.0))
        print(50, f(50.0))
        print(100, f(100.0))
        print(200, f(200.0))
    """
    norm = Normalize(vmin=float(min_val), vmax=float(max_val))
    # sm   = ScalarMappable(norm=norm, cmap=map_name)
    # sm.set_clim(vmin=min_val, vmax=max_val)
    # cmap = get_cmap()
    cmap = get_cmap(name=map_name, lut=levels)
    cmap.set_bad(bad)
    cmap.set_over(over)
    cmap.set_under(under)
    f = lambda x: cmap(norm(x))
    return norm, cmap, f



class BrowserUI(flx.PyWidget):
    genome_name      = flx.StringProp("", settable=True, doc="current genome")
    bin_width        = flx.IntProp(   -1, settable=True, doc="current bin width")
    metric           = flx.StringProp("", settable=True, doc="current metric")
    # chromosome_name  = flx.StringProp("", settable=True, doc="current chromosome")

    _combo_genome_names_text     = "Select genome"
    _combo_bin_widths_text       = "Select bin width"
    _combo_metrics_text          = "Select distance metric"

    _is_loaded = False

    def init(self, parent, height):
        super().init()

        self._is_loaded = False

        self.combo_genome_names     = ui.ComboBox(text=self._combo_genome_names_text    , editable=False, options=[], flex=1, parent=parent, css_class="combo_sel")
        self.combo_bin_widths       = ui.ComboBox(text=self._combo_bin_widths_text      , editable=False, options=[], flex=1, parent=parent, css_class="combo_sel")
        self.combo_metrics          = ui.ComboBox(text=self._combo_metrics_text         , editable=False, options=[], flex=1, parent=parent, css_class="combo_sel")

        self._is_loaded = True

        setHidden(self.combo_bin_widths      , True)
        setHidden(self.combo_metrics         , True)

    # def _render_dom(self):
    #     # Use this to determine the content. This method may return a
    #     # string, a list of virtual nodes, or a single virtual node
    #     # (which must match the type produced in _create_dom()).
    #     return [flx.create_element('span', {},
    #                 'Hello', flx.create_element('b', {}, self.name), '! '),
    #             flx.create_element('span', {},
    #                 'I happen to know that your age is %i.' % self.age),
    #             flx.create_element('br'),
    #             flx.create_element('button', {'onclick': self.increase_age},
    #                 'Next year ...')
    #             ]

    @flx.action
    def set_genome_name(self, genome_name: str):
        if self._is_loaded and genome_name != "":
            print("BrowserUI.set_genome_name", genome_name)
            self._mutate_genome_name(genome_name)
            self.reset_bin_widths()
            self.update_bin_widths()
            setHidden(self.combo_bin_widths, False)

    @flx.action
    def set_bin_width(self, bin_width: int):
        if self._is_loaded and bin_width != -1 and bin_width != "":
            print("BrowserUI.set_bin_width", bin_width)
            self._mutate_bin_width(bin_width)
            self.reset_metrics()
            self.update_metrics()
            setHidden(self.combo_metrics, False)

    @flx.action
    def set_metric(self, metric: str):
        if self._is_loaded and metric != "":
            print("BrowserUI.set_metric", metric)
            self._mutate_metric(metric)

            if (
                self._is_loaded            and
                self.genome_name     != "" and
                self.bin_width       != -1 and self.bin_width != "" and
                self.metric          != ""
            ):
                print("BrowserUI.set_metric", self.genome_name, self.bin_width, self.metric)

                self.root.genomes.update_genome(self.genome_name, self.bin_width, self.metric)

    @flx.action
    def reset_genome_names(self):
        setHidden(self.combo_genome_names, True)
        self.set_genome_names([])
        self.set_genome_names("")
        self.combo_genome_names.set_text(self._combo_genome_names_text)
        self.reset_bin_widths()

    @flx.action
    def reset_bin_widths(self):
        setHidden(self.combo_bin_widths, True)
        self.set_bin_widths([])
        self.set_bin_width(-1)
        self.combo_bin_widths.set_text(self._combo_bin_widths_text)
        self.reset_metrics()

    @flx.action
    def reset_metrics(self):
        setHidden(self.combo_metrics, True)
        self.set_metrics([])
        self.set_metric("")
        self.combo_metrics.set_text(self._combo_metrics_text)


    @flx.action
    def update_genome_names(self):
        print("BrowserUI.update_genome_names")
        val = self.root.genomes.genomes()
        self.set_genome_names(val)
        self.combo_genome_names.set_selected_index(0)

    @flx.action
    def update_bin_widths(self):
        print("BrowserUI.update_bin_widths")
        val = self.root.genomes.bin_widths(self.genome_name)
        self.set_bin_widths(val)
        self.combo_bin_widths.set_selected_index(0)

    @flx.action
    def update_metrics(self):
        print("BrowserUI.update_metrics")
        val = self.root.genomes.metrics(self.genome_name, self.bin_width)
        self.set_metrics(val)
        self.combo_metrics.set_selected_index(0)



    @flx.action
    def set_genome_names(self, val):
        print("BrowserUI.set_genome_names")
        self.combo_genome_names.set_options(val)

    @flx.action
    def set_bin_widths(self, val):
        print("BrowserUI.set_bin_widths")
        self.combo_bin_widths.set_options(val)

    @flx.action
    def set_metrics(self, val):
        print("BrowserUI.set_metrics")
        self.combo_metrics.set_options(val)



    #https://flexx.readthedocs.io/en/v0.8.0/ui/dropdown.html?highlight=dropdown
    #https://flexx.readthedocs.io/en/v0.8.0/guide/reactions.html?highlight=reaction
    @flx.reaction("combo_genome_names.selected_key")
    def reaction_combo_genome_names(self, ev):
        self.set_genome_name(self.combo_genome_names.selected_key)

    @flx.reaction("combo_bin_widths.selected_key")
    def reaction_combo_bin_widths(self, ev):
        self.set_bin_width(self.combo_bin_widths.selected_key)

    @flx.reaction("combo_metrics.selected_key")
    def reaction_combo_metrics(self, ev):
        self.set_metric(self.combo_metrics.selected_key)


    @flx.action
    def update_genome(self):
        print("BrowserUI.update_genome")
        setHidden(self.combo_chromosome_names, False)



class Graph(flx.CanvasWidget):
    CSS = """
    .flx-Drawing {
        background: #FFF;
        border: 1px solid #000;
    }
    .flx-Graph {
        width: 100% !important;
        height: 100% !important;
        overflow: scroll;
    }
    """

    def init(self, parent):
        super().init()
        self.set_capture_wheel(False)
        # self.apply_style('overflow: scroll;')  # enable scrolling
        self.ctx = self.node.getContext('2d')

    @flx.action
    def resize(self, width: int, height: int):
        print(f"Graph :: setting size :: width {width} height {height} outernode {self.outernode} node {self.node}")

        for node in [self.node]:
        # for node in [self.outernode, self.node]:
            if node:
                node.setAttribute("width" , f"{width}px")
                node.setAttribute("height", f"{height}px")
                node.style["width"]       = f"{width}px"
                node.style["height"]      = f"{height}px"
                node.style["min-width"]   = f"{width}px"
                node.style["min-height"]  = f"{height}px"
                node.style["max-width"]   = f"{width}px"
                node.style["max-height"]  = f"{height}px"

        if hasattr(self, "ctx"):
            print("Graph :: resizing canvas", width, height)
            self.ctx.width  = width
            self.ctx.height = height

    @flx.action
    def set_points(self, coord_data):
        print(f"Graph.set_points :: coord_data")

        matrix_min      = coord_data["min"            ]
        matrix_max      = coord_data["max"            ]
        img_shape       = coord_data["shape"          ]
        imgd            = coord_data["data"           ]
        colorlist       = coord_data["colorlist"      ]
        samplenames     = coord_data["samplenames"    ]
        binwidth        = coord_data["binwidth"       ]
        binnames        = coord_data["binnames"       ]
        chromosome_name = coord_data["chromosome_name"]
        vcf_name        = coord_data["vcf_name"       ]
        bin_width       = coord_data["bin_width"      ]
        metric          = coord_data["metric"         ]
        reference_name  = coord_data["sample_name"    ]
        colorrange      = coord_data["colorrange"     ]
        min_val         = coord_data["min_val"        ]
        max_val         = coord_data["max_val"        ]
        snps            = coord_data["snps"           ]
        first_position  = coord_data["first_position" ]
        last_position   = coord_data["last_position"  ]
        bin_colorrange  = coord_data["bin_colorrange" ]
        bin_snps_color  = coord_data["bin_snps_color" ]
        bin_snps        = coord_data["bin_snps"       ]
        bin_snps_min    = coord_data["bin_snps_min"   ]
        bin_snps_max    = coord_data["bin_snps_max"   ]

        num_bins        = img_shape[0]
        num_samples     = img_shape[1]
        # rgba       = img_shape[2]

        print(f"Graph.set_points :: coord_data :: matrix_min             {matrix_min}")
        print(f"Graph.set_points :: coord_data :: matrix_max             {matrix_max}")
        print(f"Graph.set_points :: coord_data :: img_shape              {img_shape}")
        print(f"Graph.set_points :: coord_data :: num_bins               {num_bins}")
        print(f"Graph.set_points :: coord_data :: num_samples            {num_samples}")
        print(f"Graph.set_points :: coord_data :: binwidth               {binwidth}")
        print(f"Graph.set_points :: coord_data :: binnames               {len(binnames)}")
        # print(f"Graph.set_points :: coord_data :: rgba       {rgba}")
        # print(f"Graph.set_points :: coord_data :: imgd       {imgd}")


        print(f"Graph.set_points :: clearing")
        self.reset()

        ctx                    = self.ctx
        
        min_display_font       = 10
        font_size              = 12
        font_height            = font_size
        font_width             = font_size * 0.48 #https://www.lifewire.com/aspect-ratio-table-common-fonts-3467385
        header_lines           = 6 + 1

        print(f"Graph.set_points :: coord_data :: font_size              {font_size}")
        print(f"Graph.set_points :: coord_data :: font_height            {font_height}")
        print(f"Graph.set_points :: coord_data :: font_width             {font_width}")

        header_height          = header_lines * font_height + font_height
        max_bin_name_length    = max([len(b) for  b in binnames])
        max_bin_name_size      = max_bin_name_length    * font_width
        max_bin_name_offset    = max_bin_name_size            + font_width + header_height
        max_bin_name_offset_h  = max_bin_name_size      * 1.2 + font_width + header_height
        print(f"Graph.set_points :: coord_data :: max_bin_name_length    {max_bin_name_length}")
        print(f"Graph.set_points :: coord_data :: max_bin_name_size      {max_bin_name_size}")
        print(f"Graph.set_points :: coord_data :: max_bin_name_offset    {max_bin_name_offset}")

        max_sample_name_length = max([len(x) for x in samplenames])
        max_sample_name_size   = max_sample_name_length * font_width
        print(f"Graph.set_points :: coord_data :: max_sample_name_length {max_sample_name_length}")
        print(f"Graph.set_points :: coord_data :: max_sample_name_size   {max_sample_name_size}")

        start_y                = (max_bin_name_offset  * (font_size >= min_display_font)) +   font_width
        start_x                = (max_sample_name_size * (font_size >= min_display_font)) + 4*font_width
        print(f"Graph.set_points :: coord_data :: start_y                {start_y}")
        print(f"Graph.set_points :: coord_data :: start_x                {start_x}")

        height_eff             = start_y + num_samples * font_height + font_height
        width_eff              = start_x + num_bins    * font_height + font_height
        print(f"Graph.set_points :: coord_data :: height_eff             {height_eff}")
        print(f"Graph.set_points :: coord_data :: width_eff              {width_eff}")

        print(f"Graph.set_points :: resize :: width {width_eff} height {height_eff}")
        self.resize(width_eff, height_eff)

        ##
        ## Write Border
        ##
        if True:
            print(f"Graph.set_points :: Writing Border")
            ctx.strokeStyle = 'black'
            # ctx.fillStyle   = "rgba(255, 255, 255, 0.0)"
            ctx.lineWidth   = 1
            ctx.beginPath()
            ctx.rect(
                1,
                1,
                width_eff  - 2,
                height_eff - 2
            )
            ctx.stroke()
            ctx.lineWidth   = 0

        ##
        ## Write title
        ##
        if font_size >= min_display_font:
            print(f"Graph.set_points :: Writing Title")

            ctx.font        = f"{font_size}px Courier New"
            ctx.strokeStyle = 'black'
            ctx.fillStyle   = 'black'
            ctx.save()
            ctx.translate(start_x + font_height, max_bin_name_offset_h)
            ctx.rotate(-0.5*math.pi)
            for bin_pos, bin_name in enumerate(binnames):
                text_x = 0
                text_y = (bin_pos * font_height)
                ctx.fillText(bin_name, text_x, text_y)
            ctx.restore()


            # ctx.font        = f"{font_size}px Arial"
            # ctx.fillStyle   = "black"
            # ctx.strokeStyle = "black"
            # ctx.textAlign   = "left"
            for sample_pos, sample_name in enumerate(samplenames):
                text_x = font_width
                text_y = (sample_pos * font_height) + start_y + font_height
                ctx.fillText(sample_name, text_x, text_y)

            ##
            ## Header
            ##   header_lines
            if True:
                print(f"Graph.set_points :: Writing Header")
                ctx.save()

                title_max_len    = max([len(v) for v in [chromosome_name, vcf_name, bin_width, metric, reference_name, first_position, last_position, snps]])
                title_fmt        = "{:"+str(title_max_len)+"}"

                text_x  = font_width
                text_y  = 1 * font_height
                ctx.fillText(("VCF           : "+title_fmt).format(vcf_name         ), text_x, text_y)

                text_x  = font_width
                text_y  = 2 * font_height
                ctx.fillText(("Chromosome    : "+title_fmt).format(chromosome_name  ), text_x, text_y)
                text_x += (title_max_len + 16) * font_width
                text_y  = 2 * font_height
                ctx.fillText(("Bin Width     : "+title_fmt).format(bin_width        ), text_x, text_y)

                text_x  = font_width
                text_y  = 3 * font_height
                ctx.fillText(("Metric        : "+title_fmt).format(metric           ), text_x, text_y)
                text_x += (title_max_len + 16) * font_width
                text_y  = 3 * font_height
                ctx.fillText(("Reference     : "+title_fmt).format(reference_name   ), text_x, text_y)

                text_x  = font_width
                text_y  = 4 * font_height
                ctx.fillText(("First Position: "+title_fmt).format(first_position   ), text_x, text_y)
                text_x += (title_max_len + 16) * font_width
                text_y  = 4 * font_height
                ctx.fillText(("Last Position : "+title_fmt).format(last_position    ), text_x, text_y)

                text_x  = font_width
                text_y  = 5 * font_height
                ctx.fillText(("# SNPs        : "+title_fmt).format(snps             ), text_x, text_y)

                ctx.restore()

            ##
            ## Color list
            ##
            if True:
                print(f"Graph.set_points :: Writing Color List")
                blocksize  = 1
                block_y    = font_height * 6
                block_x_ts = font_width
                block_x_bs = block_x_ts + 6*font_width
                block_x_te = block_x_bs + (len(colorrange) * blocksize) + (1*font_width)
                ctx.save()
                ctx.fillText("{:3.2f}".format(min_val), block_x_ts, block_y)
                ctx.fillText("{:3.2f}".format(max_val), block_x_te, block_y)
                for color_num, color_pos in enumerate(colorrange):
                    color           = colorlist[color_pos]
                    ctx.strokeStyle = color
                    ctx.fillStyle   = color
                    ctx.beginPath()
                    ctx.fillRect(
                        block_x_bs  + (color_num * blocksize),
                        block_y - font_height*.7,
                        blocksize,
                        font_height
                    )
                    ctx.stroke()
                ctx.restore()

        # bin_colorrange  = coord_data["bin_colorrange" ]
        # bin_snps_color  = coord_data["bin_snps_color" ]
        # bin_snps        = coord_data["bin_snps"       ]
        # bin_snps_min    = coord_data["bin_snps_min"   ]
        # bin_snps_max    = coord_data["bin_snps_max"   ]

        ##
        ## Create Graph
        ##
        if True:
            print(f"Graph.set_points :: creating graph")
            print(f"Graph.set_points :: setting")
            print(f"Graph.set_points :: setting :: imgd.length {imgd.length}")
            print(f"Graph.set_points :: setting :: imgd[0].length {imgd[0].length}")
            for x in range(num_bins):
                col = imgd[x]
                for y in range(num_samples):
                    colorpos        = col[y]
                    color           = colorlist[colorpos]
                    ctx.strokeStyle = color
                    ctx.fillStyle   = color
                    ctx.beginPath()
                    ctx.fillRect(
                        start_x + x * font_height,
                        start_y + y * font_height,
                        font_height,
                        font_height
                    )
                    ctx.stroke()
            print(f"Graph :: done")
        
        return
        if False:
            tray_rows    = coord_data["tray_rows"]
            tray_columns = coord_data["tray_columns"]
            image_width  = coord_data["image_width"]
            image_height = coord_data["image_height"]
            marked_cells = coord_data["marked_cells"]
            positions    = coord_data["positions"]
            shadows      = coord_data["shadows"]

            # print("Drawing.set_points", width, height)

            # self.label.set_text(t.format(ev.pos[0], ev.pos[1], ev.buttons))
            if any([v == -1 for v in [image_width, image_height, self.scale_factor]]):
                print("Drawing.set_points - not set")
                return

            if positions is None:
                print("Drawing.set_points - no positions")
                return

            scale_factor    = self.scale_factor/100
            border          = 1
            lineWidth       = 3
            gapBetweenCells = 3

            width_eff       = image_width  * scale_factor
            height_eff      = image_height * scale_factor
            ctx             = self.ctx

            self.resize(width_eff + 4*border, height_eff + 4*border)

            # print(f"Drawing.set_points - clearing - width {self.width} height {self.height} scale_factor {scale_factor}")
            # print(f"Drawing.set_points - clearing - {self.width  * scale_factor} {self.height * scale_factor}")

            ctx.clearRect(
                0,
                0,
                width_eff  + 4*border,
                height_eff + 4*border
            )

            if self.show_grid:
                ctx.lineWidth   = lineWidth
                cell_width      = width_eff  // tray_columns
                cell_height     = height_eff // tray_rows
                for y_num, square_y_pos in enumerate(range(0, height_eff, cell_height)):
                    for x_num, square_x_pos in enumerate(range(0, width_eff, cell_width)):
                        if (y_num, x_num) in marked_cells:
                            ctx.strokeStyle = "#089000"
                        else:
                            ctx.strokeStyle = "#FF0000"

                        ctx.beginPath()
                        ctx.rect(
                            square_x_pos + gapBetweenCells,
                            square_y_pos + gapBetweenCells,
                            cell_width   - gapBetweenCells - lineWidth//2,
                            cell_height  - gapBetweenCells - lineWidth//2
                        )
                        ctx.stroke()

            if self.show_shadows:
                ctx.strokeStyle = self.shadowStrokeStyle
                ctx.fillStyle   = self.shadowFillStyle

                for pos in range(len(shadows[0])):
                    y, x = shadows[0][pos], shadows[1][pos]
                    # print(f"Drawing.set_points - x {x} y {y} X {x * self.scale_factor} Y {y * self.scale_factor} scale_factor {self.scale_factor} {self.strokeStyle} {self.fillStyle}")
                    ctx.beginPath()
                    ctx.fillRect(
                        x * scale_factor + 1,
                        y * scale_factor + 1,
                        scale_factor,
                        scale_factor
                    )
                    ctx.stroke()

            ctx.strokeStyle = self.highlightStrokeStyle
            ctx.fillStyle   = self.highlightFillStyle

            for pos in range(len(positions[0])):
                y, x = positions[0][pos], positions[1][pos]
                # print(f"Drawing.set_points - x {x} y {y} X {x * self.scale_factor} Y {y * self.scale_factor} scale_factor {self.scale_factor} {self.strokeStyle} {self.fillStyle}")
                ctx.beginPath()
                ctx.fillRect(
                    x * scale_factor + 1,
                    y * scale_factor + 1,
                    scale_factor,
                    scale_factor
                )
                ctx.stroke()

    @flx.action
    def reset(self, zero=True):
        print("Graph.reset", zero)
        if hasattr(self.ctx, 'self.width'):
            print("Graph.reset :: clearing")
            self.ctx.clearRect(
                0,
                0,
                self.ctx.width,
                self.ctx.height
            )
        if zero:
            print("Graph.reset :: zeroing")
            self.resize(0, 0)


class ChromosomeController(flx.PyWidget):
    filename                  = flx.StringProp(    "-", settable=True, doc="file name")

    vcf_name                  = flx.StringProp(    "-", settable=True, doc="vcf name")
    bin_width                 = flx.IntProp(        -1, settable=True, doc="current bin width")
    chromosome_order          = flx.IntProp(        -1, settable=True, doc="order of the chromosome")
    chromosome_name           = flx.StringProp(    "-", settable=True, doc="name of the chromosome")
    metric                    = flx.StringProp(    "-", settable=True, doc="statistics metric")

    matrix_size               = flx.IntProp(        -1, settable=True, doc="size of triangular matrix")
    bin_max                   = flx.IntProp(        -1, settable=True, doc="minimum bin value")
    bin_min                   = flx.IntProp(        -1, settable=True, doc="maximum bin value")
    bin_count                 = flx.IntProp(        -1, settable=True, doc="number of bins")

    bin_snps_min              = flx.IntProp(        -1, settable=True, doc="minimum number of snps in a bin")
    bin_snps_max              = flx.IntProp(        -1, settable=True, doc="maximum number of snps in a bin")
    bin_snps                  = flx.ListProp(       [], settable=True, doc="number of snps per bin")

    chromosome_snps           = flx.IntProp(        -1, settable=True, doc="number of snps in the chromosome")
    chromosome_first_position = flx.IntProp(        -1, settable=True, doc="first chromosome position")
    chromosome_last_position  = flx.IntProp(        -1, settable=True, doc="last chromosome position")

    sample_names              = flx.ListProp(       [], settable=True, doc="sample names")
    sample_count              = flx.IntProp(        -1, settable=True, doc="number of samples")
    sample_name               = flx.StringProp(    "-", settable=True, doc="sample name")

    min_val                   = flx.FloatProp(     0.0, settable=True, doc="minimun bin value")
    max_val                   = flx.FloatProp(     1.0, settable=True, doc="maximum bin value")
    num_vals                  = flx.IntProp(       256, settable=True, doc="number of steps")

    color_name                = flx.StringProp("Greys", settable=True, doc="color map")
    color_bad                 = flx.StringProp("pink" , settable=True, doc="color map")
    color_over                = flx.StringProp("pink" , settable=True, doc="color map")
    color_under               = flx.StringProp("pink" , settable=True, doc="color map")

    _is_loaded = False

    def init(self, parent, height):
        super().init()

        self._is_loaded = False

        self.chromosome : reader.Chromosome = None

        self.cmap = None

        colors = []
        for k,v in color_maps.items():
            for z in v:
                colors.append((z,f"{k} - {z}"))

        with flx.VBox(flex=4, parent=parent):
            with flx.HBox(flex=1):
                flx.Label(text=lambda: f"Database Filename: {os.path.basename(self.filename)}", flex=1)
                flx.Label(text=lambda: f"# Matrix Size: {self.matrix_size:9,d}"               , flex=1)
                flx.Label(text=lambda: f"First Position: {self.chromosome_first_position:9,d}", flex=1)
                flx.Label(text=lambda: f"Last Position: {self.chromosome_last_position:9,d}"  , flex=1)

            with flx.HBox(flex=1):
                flx.Label(text=lambda: f"# Bins: {self.bin_count:9,d}"                        , flex=1)
                flx.Label(text=lambda: f"# Bin Min: {self.bin_min:9,d}"                       , flex=1)
                flx.Label(text=lambda: f"# Bin Max: {self.bin_max:9,d}"                       , flex=1)

            with flx.HBox(flex=1):
                flx.Label(text=lambda: f"# SNPs: {self.chromosome_snps:9,d}"                  , flex=1)
                flx.Label(text=lambda: f"# Bin SNPs Min: {self.bin_snps_min:9,d}"             , flex=1)
                flx.Label(text=lambda: f"# Bin SNPs Max: {self.bin_snps_max:9,d}"             , flex=1)

            with flx.HBox(flex=1):
                self.combo_color_names = ui.ComboBox(text="Select color scheme", editable=False, options=colors, flex=1)
                self.combo_color_names.set_text(self.color_name)
                # flx.Label(text=lambda: f"Color name: {self.color_name}"             , flex=1)

                self.combo_sample_names = ui.ComboBox(text="Select sample name", editable=False, options=[], flex=1)
                self.combo_sample_names.set_text("Select sample name")
                # flx.Label(text=lambda: f"Sample name: {self.sample_name}"             , flex=1)

            with flx.HBox(flex=height-4) as graphbox:
                self.graph = Graph(graphbox)

        self._is_loaded = True

    @flx.action
    def reset(self):
        print("ChromosomeController.reset")
        self.chromosome = None

        self.set_filename(                 "-")
    
        self.set_vcf_name(                 "-")
        self.set_bin_width(                 -1)
        self.set_chromosome_order(          -1)
        self.set_chromosome_name(          "-")
        self.set_metric(                   "-")

        self.set_matrix_size(               -1)
        self.set_bin_max(                   -1)
        self.set_bin_min(                   -1)
        self.set_bin_count(                 -1)

        self.set_bin_snps_min(              -1)
        self.set_bin_snps_max(              -1)
        self.set_bin_snps(                  [])

        self.set_chromosome_snps(           -1)
        self.set_chromosome_first_position( -1)
        self.set_chromosome_last_position(  -1)

        self.set_sample_names(              [])
        self.set_sample_count(              -1)
        self.set_sample_name(              "-")

        self.combo_sample_names.set_text("Select sample name")
        self.combo_sample_names.set_options(self.sample_names)
        self.combo_sample_names.set_selected_index(None)
        self.graph.reset()
        print("ChromosomeController.reset :: done")

    @flx.action
    def update(self):
        print("ChromosomeController.update")

        if self.chromosome is None:
            print("ChromosomeController.update :: NO CHROMOSOME")

        else:
            print("ChromosomeController.update :: updating")
            self.set_filename(                 self.chromosome.filename)

            self.set_vcf_name(                 self.chromosome.vcf_name)
            self.set_bin_width(                self.chromosome.bin_width)
            self.set_chromosome_order(         self.chromosome.chromosome_order)
            self.set_chromosome_name(          self.chromosome.chromosome_name)
            self.set_metric(                   self.chromosome.metric)

            self.set_matrix_size(              self.chromosome.matrix_size)
            self.set_bin_max(                  self.chromosome.bin_max)
            self.set_bin_min(                  self.chromosome.bin_min)
            self.set_bin_count(                self.chromosome.bin_count)

            self.set_bin_snps_min(             self.chromosome.bin_snps_min)
            self.set_bin_snps_max(             self.chromosome.bin_snps_max)
            self.set_bin_snps(                 self.chromosome.binsnpNp.tolist())

            self.set_chromosome_snps(          self.chromosome.chromosome_snps)
            self.set_chromosome_first_position(self.chromosome.chromosome_first_position)
            self.set_chromosome_last_position( self.chromosome.chromosome_last_position)

            self.set_sample_names(             self.chromosome.sample_names)
            self.set_sample_count(             self.chromosome.sample_count)
            self.set_sample_name(              "-")

            self.combo_sample_names.set_text("Select sample name")
            self.combo_sample_names.set_options(self.sample_names)
            self.combo_sample_names.set_selected_index(0)

    @flx.action
    def set_chromosome(self, chromosome: reader.Chromosome):
        print("ChromosomeController.set_chromosome")

        self.reset()

        self.chromosome = chromosome

        if self.chromosome is not None:
            self.update()

    @flx.reaction
    def reaction_min_val(self):#, ev):
        self.cmap, _ = self.update_color()

    @flx.reaction
    def reaction_max_val(self):#, ev):
        self.cmap, _ = self.update_color()

    @flx.reaction("combo_color_names.selected_key")
    def reaction_combo_color_names(self, *args):
        print("ChromosomeController.reaction_combo_color_names", self.combo_color_names.selected_key)
        self.set_color_name(self.combo_color_names.selected_key)

    @flx.reaction("combo_sample_names.selected_key")
    def reaction_combo_sample_names(self, *args):
        print("ChromosomeController.reaction_combo_sample_names", self.combo_sample_names.selected_key)
        self.set_sample_name(self.combo_sample_names.selected_key)

    @flx.reaction
    def reaction_sample_name(self):
        print("ChromosomeController.reaction_sample_name", self.sample_name)
        self.display()

    def update_color(self, color_name=None, min_val=None, max_val=None):
        color_name = self.color_name if color_name is None else color_name
        min_val    = self.min_val    if min_val    is None else min_val
        max_val    = self.max_val    if max_val    is None else max_val
        val_ran    = (max_val-min_val)/self.num_vals
        _, _, cmap = get_color_maps(
            color_name,
            min_val = min_val,
            max_val = max_val,
            levels  = self.num_vals,
            bad     = self.color_bad,
            over    = self.color_over,
            under   = self.color_under
        )
        arange = np.arange(min_val, max_val, val_ran)
        acolor = cmap(arange)
        # print("min_val ", min_val)
        # print("max_val ", max_val)
        # print("val_ran ", val_ran)
        # print("num_vals", self.num_vals)
        # print("arange  ", arange)
        # print("acolor  ", acolor)
        return cmap, acolor

    def display(self):
        print(f"ChromosomeController.display :: _is_loaded {self._is_loaded} sample_name {self.sample_name} color_name {self.color_name}")
        if not (
            self.chromosome      is not None and
            self._is_loaded                  and
            self.chromosome_name is not None and self.chromosome_name != "-" and self.chromosome_name != "" and
            self.sample_name     is not None and self.sample_name     != "-" and self.sample_name     != "" and
            self.color_name      is not None and self.color_name      != "-" and self.color_name      != ""
        ):
            print(f"ChromosomeController.display :: INCOMPLETE")
        else:
            print(f"ChromosomeController.display :: displaying")
            matrix     = self.matrix_sample(self.sample_name)
            # print(f"ChromosomeController.display :: matrix", matrix)
            # print(f"ChromosomeController.display :: np.nanmin(matrix)", np.nanmin(matrix))
            # print(f"ChromosomeController.display :: np.nanmax(matrix)", np.nanmax(matrix))
            matrix     = matrix.astype(np.float)
            matrix[matrix == np.nanmin(matrix)] = None
            matrix[matrix == np.nanmax(matrix)] = None
            matrix_min = float(np.nanmin(matrix))
            matrix_max = float(np.nanmax(matrix))

            bin_snps_min = min(self.bin_snps)
            bin_snps_max = max(self.bin_snps)
            # bin_snps_pro = [(b-bin_snps_min)/(bin_snps_max-bin_snps_min) for b in self.bin_snps]
            # print("self.bin_snps", self.bin_snps)
            # print("bin_snps_min ", bin_snps_min )
            # print("bin_snps_max ", bin_snps_max )
            # print("bin_snps_pro ", bin_snps_pro )
            # print(list(zip(self.bin_snps, bin_snps_pro)))
            bin_cmap, bin_acolor = self.update_color(min_val=bin_snps_min, max_val=bin_snps_max)
            bin_cmapg  = bin_cmap(self.bin_snps)
            bin_cmapg  = bin_cmapg.tolist()
            bin_acolor = bin_acolor.tolist()
            # print("bin_cmapg", bin_cmapg)

            # print(f"ChromosomeController.display :: matrix = {matrix}")
            print(f"ChromosomeController.display :: matrix max = {matrix_max}")
            print(f"ChromosomeController.display :: matrix min = {matrix_min}")
            self.set_min_val(matrix_min)
            self.set_max_val(matrix_max)

            self.cmap, acolor = self.update_color(min_val=matrix_min, max_val=matrix_max)
            img               = self.cmap(matrix)

            # print(f"ChromosomeController.display :: img = {img}")
            print(f"ChromosomeController.display :: img shape  = {img.shape}")
            print(f"ChromosomeController.display :: img type   = {type(img)}")
            print(f"ChromosomeController.display :: bin_count  = {self.bin_count}")
            # print(f"ChromosomeController.display :: bin_count  = {self.bin_n}")
            imgd      = img.tolist()
            acolord   = acolor.tolist()
            colorlist = {}

            def col_to_color(col):
                for row_num in range(len(col)):
                    cell     = col[row_num]
                    rgba     = [int(v*255) for v in cell]
                    r,g,b,a  = rgba
                    hexa     = f"#{r:02x}{g:02x}{b:02x}{a:02x}"
                    colorpos = colorlist.get(hexa, -1)
                    if colorpos == -1:
                        colorpos = len(colorlist)
                        colorlist[hexa] = colorpos
                    col[row_num] = colorpos

            for col_num in range(len(imgd)):
                col = imgd[col_num]
                col_to_color(col)
            
            col_to_color(acolord)
            # print("acolord", acolord, type(acolord))
            col_to_color(bin_cmapg)
            # print("bin_cmapg", bin_cmapg, type(bin_cmapg))
            col_to_color(bin_acolor)
            # print("bin_acolor", bin_acolor, type(bin_acolor))


            colorlist = sorted(list(colorlist.keys()), key=lambda x: colorlist[x])
            # print(f"ChromosomeController.display :: imgd = {imgd}")
            # bin_names = [f"{b+1:7,d} - {b*self.bin_width:12,d}-{(b+1)*self.bin_width:12,d}" for b in range(self.bin_count)]
            bin_names = self.format_range()
            self.graph.set_points({
                "min"            : matrix_min,
                "max"            : matrix_max,
                "shape"          : img.shape,
                "data"           : imgd,
                "colorlist"      : colorlist,
                "samplenames"    : self.sample_names,
                "binwidth"       : self.bin_width,
                "binnames"       : bin_names,
                "chromosome_name": self.chromosome_name,
                "vcf_name"       : os.path.basename(self.vcf_name),
                "bin_width"      : f"{self.bin_width:12,d}",
                "metric"         : self.metric,
                "sample_name"    : self.sample_name,
                "colorrange"     : acolord,
                "min_val"        : matrix_min,
                "max_val"        : matrix_max,
                "snps"           : f"{self.chromosome_snps:12,d}",
                "first_position" : f"{self.chromosome_first_position:12,d}",
                "last_position"  : f"{self.chromosome_last_position:12,d}",
                "bin_colorrange" : bin_acolor,
                "bin_snps_color" : bin_cmapg,
                "bin_snps"       : self.bin_snps,
                "bin_snps_min"   : bin_snps_min,
                "bin_snps_max"   : bin_snps_max

            })

    def format_range(self):
        def human_format(num, suffixes=[['bp',0], ['Kbp',0], ['Mbp',2], ['Gbp',3], ['Tbp',4], ['Pbp',5]]):
            # m = sum([abs(num/1000.0**x) >= 1 for x in range(1, len(suffixes))])
            m = int(math.log10(num) // 3) if num > 0 else 0
            return f'{num/1000.0**m:.{suffixes[m][1]}f}{suffixes[m][0]}'

        res = [None] * self.bin_count
        for b in range(self.bin_count):
            i,s,e = b, b*self.bin_width, (b+1)*self.bin_width
            s     = human_format(s)
            e     = human_format(e)
            res[b] = f"{i+1:6,d}:{s}-{e}"
            # bin_names = [f"{b+1:7,d} - {b*self.bin_width:12,d}-{(b+1)*self.bin_width:12,d}" for b in range(self.bin_count)]
        return res

    def matrix_sample(self, sample_name) -> np.ndarray:
        return self.chromosome.matrix_sample(sample_name, metric=self.metric)

    def matrix_bin(self, binNum: int) -> np.ndarray:
        return self.chromosome.matrix_bin(binNum)

    def matrix_bin_square(self, binNum: int) -> np.ndarray:
        return self.chromosome.matrix_bin_square(binNum)

    def matrix_bin_sample(self, binNum: int, sample_name: str) -> np.ndarray:
        return self.chromosome.matrix_bin_sample(binNum, sample_name)

    def matrix_bin_dist(self, binNum: int) -> np.ndarray:
        return self.chromosome.matrix_bin_dist(binNum, metric=self.metric, dtype=self.matrix_dtype)

    def matrix_bin_dist_square(self, binNum: int) -> np.ndarray:
        return self.chromosome.matrix_bin_dist_square(binNum, self.metric)

    def matrix_bin_dist_sample(self, binNum: int, sample_name: str) -> np.ndarray:
        return self.chromosome.matrix_bin_dist_sample(binNum, sample_name, metric=self.metric)

    def matrix_bin_dist_sample_square(self, binNum: int, sample_name: str) -> np.ndarray:
        return self.chromosome.matrix_bin_dist_sample_square(binNum, sample_name, metric=self.metric)



class GenomeController(flx.PyWidget):
    genome_name      = flx.StringProp("-", settable=True, doc="genome name")
    vcf_name         = flx.StringProp("-", settable=True, doc="vcf name")
    bin_width        = flx.IntProp(    -1, settable=True, doc="current bin width")
    metric           = flx.StringProp("-", settable=True, doc="statistics metric")
    sample_names     = flx.ListProp([]   , settable=True, doc="sample names")
    sample_count     = flx.IntProp(    -1, settable=True, doc="number of samples")
    chromosome_names = flx.ListProp([]   , settable=True, doc="chromosome names")
    chromosome_count = flx.IntProp(    -1, settable=True, doc="number of chromosomes")
    genome_bins      = flx.IntProp(    -1, settable=True, doc="number of bins")
    genome_snps      = flx.IntProp(    -1, settable=True, doc="number of snps")
    filename         = flx.StringProp("-", settable=True, doc="filename")

    chromosome_name  = flx.StringProp("", settable=True, doc="current chromosome")

    _combo_chromosome_names_text = "Select chromosome"
    _is_loaded                   = False

    def init(self, parent, height):
        super().init()

        self._is_loaded = False

        self.genome     : reader.Genome        = None

        with flx.VBox(flex=height, parent=parent):
            with flx.HBox(flex=4):
                flx.Label(text=lambda: f"VCF Name: {os.path.basename(self.vcf_name)}"             , flex=1)
                flx.Label(text=lambda: f"Database Filename: {os.path.basename(self.filename)}"    , flex=1)
                flx.Label(text=lambda: f"Bin Width: {self.bin_width:9,d}"           , flex=1)
                flx.Label(text=lambda: f"Metric: {self.metric}"                     , flex=1)

            with flx.HBox(flex=4):
                flx.Label(text=lambda: f"# Samples: {self.sample_count:9,d}"        , flex=1)
                flx.Label(text=lambda: f"# Bins: {self.genome_bins:9,d}"            , flex=1)
                flx.Label(text=lambda: f"# SNPs: {self.genome_snps:9,d}"            , flex=1)

            with flx.HBox(flex=4):
                flx.Label(text=lambda: f"# Chromosomes: {self.chromosome_count:9,d}", flex=1)
                self.combo_chromosome_names = ui.ComboBox(text=self._combo_chromosome_names_text, editable=False, options=[], flex=3, css_class="combo_sel")

        self._is_loaded = True

        setHidden(self.combo_chromosome_names, True)

    @flx.action
    def reset(self):
        print("GenomeController.reset")
        self.genome = None
        self.set_genome_name(     "-")
        self.set_vcf_name(        "-")
        self.set_bin_width(       -1)
        self.set_metric(          "-")
        self.set_sample_names(    [])
        self.set_sample_count(    -1)
        self.set_chromosome_names([])
        self.set_chromosome_count(-1)
        self.set_genome_bins(     -1)
        self.set_genome_snps(     -1)
        self.set_filename(        "-")

        print("GenomeController.reset :: reseting chromosome name")
        self.root.chromosome.reset()
        # setHidden(self.combo_chromosome_names, True)
        self.combo_chromosome_names.set_text(self._combo_chromosome_names_text)
        self.combo_chromosome_names.set_options([])
        self.set_chromosome_name("-")
        print("GenomeController.reset :: reseting chromosome name :: done")

    @flx.action
    def update(self):
        print("GenomeController.update")

        if self.genome is not None and self.genome_name != "-":
            print("GenomeController.update :: updating :: chromosome_names", self.genome.chromosome_names)
            self.set_vcf_name(        self.genome.vcf_name)
            self.set_bin_width(       self.genome.bin_width)
            self.set_metric(          self.genome.metric)
            self.set_sample_names(    self.genome.sample_names)
            self.set_sample_count(    self.genome.sample_count)
            self.set_chromosome_names(self.genome.chromosome_names)
            self.set_chromosome_count(self.genome.chromosome_count)
            self.set_genome_bins(     self.genome.genome_bins)
            self.set_genome_snps(     self.genome.genome_snps)
            self.set_filename(        self.genome.filename)

            print("GenomeController.update :: updating :: showing chromosome names")
            self.combo_chromosome_names.set_text(self._combo_chromosome_names_text)
            self.combo_chromosome_names.set_options(self.chromosome_names)
            setHidden(self.combo_chromosome_names, False)
            print("GenomeController.update :: updating :: selecting chromosome name")
            self.combo_chromosome_names.set_selected_index(0)
            print("GenomeController.update :: updating :: done")
            # setHidden(self.combo_chromosome_names, False)

    @flx.action
    def set_genome(self, genome_name: str, genome: reader.Genome):
        if self._is_loaded:
            print("GenomeController.set_genome :: genome_name", genome_name)
            
            self.reset()

            self.genome = genome

            if self.genome is not None:
                print("GenomeController.set_genome :: valid")
                self.set_genome_name(     genome_name)
                self.update()
                print("GenomeController.set_genome :: valid :: done")

    @flx.action
    def set_chromosome_name(self, chromosome_name: str):
        if self._is_loaded:
            print("GenomeController.set_chromosome_name", chromosome_name)

            self._mutate_chromosome_name(chromosome_name)

            if (
                self.genome_name     != "" and self.genome_name     != "-" and
                self.bin_width       != -1 and self.bin_width       != "-" and
                self.metric          != "" and self.metric          != "-" and
                self.chromosome_name != "" and self.chromosome_name != "-"
            ):
                print("GenomeController.set_chromosome_name", self.genome_name, self.bin_width, self.metric, self.chromosome_name)

                self.root.genomes.update_chromosome(self.genome_name, self.bin_width, self.metric, self.chromosome_name)
            else:
                print("GenomeController.set_chromosome_name", chromosome_name, "INCOMPLETE")

    @flx.reaction("combo_chromosome_names.selected_key")
    def reaction_combo_chromosome_names(self, *args):
        print("GenomeController.reaction_combo_chromosome_names", self.chromosome_names)
        self.set_chromosome_name(self.combo_chromosome_names.selected_key)



class GenomesController(flx.PyComponent):
    def init(self):
        super().init()

        self._genomes : reader.Genomes   = None

    def reset(self):
        self.root.genome.reset()

    def set_genomes(self, genomes: reader.Genomes):
        print(f"GenomesController.set_genomes")
        self._genomes = genomes
        self.reset()

    def update_genome(self, genome_name: str, bin_width: int, metric: str):
        print(f"GenomesController.update_genome genome_name {genome_name} bin_width {bin_width} metric {metric}")
        genome = self.load_genome(genome_name, bin_width, metric)
        self.root.genome.set_genome(genome_name, genome)

    def update_chromosome(self, genome_name: str, bin_width: int, metric: str, chromosome_name: str):
        print(f"GenomesController.update_chromosome genome_name {genome_name} bin_width {bin_width} metric {metric} chromosome_name {chromosome_name}")
        chromosome = self.load_chromosome(genome_name, bin_width, metric, chromosome_name)
        self.root.chromosome.set_chromosome(chromosome)

    def genomes(self):
        return self._genomes.genomes

    def genome_info(self, genome_name: str):
        return self._genomes.genome_info(genome_name)

    def bin_widths(self, genome_name: str) -> typing.List[str]:
        return self._genomes.bin_widths(genome_name)

    def bin_width_info(self, genome_name: str, bin_width: int):
        return self._genomes.bin_width_info(genome_name, bin_width)

    def metrics(self, genome_name: str, bin_width: int) -> typing.List[str]:
        return self._genomes.metrics(genome_name, bin_width)

    def metric_info(self, genome_name: str, bin_width: int, metric: str) -> typing.List[str]:
        return self._genomes.metric_info(genome_name, bin_width, metric)

    def chromosome_names(self, genome_name: str, bin_width: int, metric: str) -> typing.List[typing.Tuple[int, str]]:
        val = self._genomes.chromosome_names(genome_name, bin_width, metric)
        val = [v[1] for v in val]
        return val

    def update(self, verbose=False):
        self._genomes.update(verbose=verbose)

    def load_genome(self, genome_name: str, bin_width: int, metric: str) -> reader.Genome:
        return self._genomes.load_genome(genome_name, bin_width, metric)

    def load_chromosome(self, genome_name: str, bin_width: int, metric: str, chromosome_name: str) -> reader.Chromosome:
        return self._genomes.load_chromosome(genome_name, bin_width, metric, chromosome_name)



class MainController(flx.PyComponent):
    CSS = """
    .flx-ComboBox {
        background: #9d9 !important;
    }
    .flx-LineEdit {
        border: 2px solid #9d9;
    }
    .flx-ComboBox.hidden {
        visibility: hidden;
    }
    //.flx-ComboBox.combo_sel {
    //    width: 250px !important;
    //    height: 50px !important;
    //}
    """

    genomes_inst : reader.Genomes = None

    # https://flexx.readthedocs.io/en/stable/examples/send_data_src.html

    def init(self):
        super().init()

        with flx.VBox(flex=1):
            with flx.HFix(flex=1) as self.header:
                self.BrowserUI  = BrowserUI(self.header, 1)
            with flx.HFix(flex=3) as self.genomeBox:
                self.genome     = GenomeController(self.genomeBox, 3)
            with flx.HFix(flex=56) as self.chromosomeBox:
                self.chromosome = ChromosomeController(self.chromosomeBox, 37)

        self.genomes    = GenomesController()

        self.genomes.set_genomes(MainController.genomes_inst)
        self.genomes.update(verbose=True)
        self.BrowserUI.update_genome_names()



def main():
    folder_name = sys.argv[1]

    MainController.genomes_inst = reader.Genomes(folder_name)

    # https://flexx.readthedocs.io/en/stable/guide/running.html
    # https://flexx.readthedocs.io/en/stable/guide/reactions.html

    flexx.config.hostname           = '0.0.0.0'
    flexx.config.port               = 5000
    flexx.config.log_level          = "DEBUG" if DEBUG else "INFO"
    flexx.config.tornado_debug      = DEBUG
    flexx.config.ws_timeout         = 20
    flexx.config.browser_stacktrace = True
    flexx.config.cookie_secret      = "0123456789"

    app = flx.App(MainController)
    app.serve('')  # Serve at http://domain.com/
    flx.start()  # mainloop will exit when the app is closed

if __name__ == "__main__":
    main()
