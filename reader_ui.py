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
from pscript import RawJS

import reader

DEBUG       = True


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
        .graph {
            top: 35px;
        }
    """

    def init(self, parent):
        super().init()

        self.set_capture_wheel(False)
        # self.apply_style('overflow: scroll;')  # enable scrolling
        self.ctx = self.node.getContext('2d')
        self.set_css_class("graph")

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



class Forms(flx.Widget):
    # https://flexx.readthedocs.io/en/stable/guide/widget_basics.html
    CSS = """
        .selector:hover {
            height: 33%;
        }
        .selector {
            position:absolute;
            top: 5px;
            left: 5px;
            background-color: darkseagreen;
            width: 99%;
            height: 30px;
            z-index: 1000;
        }

        .form_sel {
            background-color: lightgreen;
            border: 1px solid #000;
            width: 15%;
            font-size: 12px;
            height: 25px;
        }

        .btn_download {
            background-color: lightgreen;
            width: 10%;
            height: 26px;
        }

        .info {
            font-size: small;
        }

        .info_selection_sp{
            width: 15%;
            display: inline-block;
        }
        .info_selection_he{
            font-weight: bold;
        }
        .info_selection_va{
        }

        .info_genome_sp{
            width: 32%;
            display: inline-block;
        }
        .info_genome_he{
            font-weight: bold;
        }
        .info_genome_va{
        }

        .info_chromosome_sp{
            width: 32%;
            display: inline-block;
        }
        .info_chromosome_he{
            font-weight: bold;
        }
        .info_chromosome_va{
        }
    """

    genomes          = flx.ListProp(  [] , settable=True, doc="List of genomes")
    bin_widths       = flx.ListProp(  [] , settable=True, doc="List of bin width")
    metrics          = flx.ListProp(  [] , settable=True, doc="List of metrics")
    chromosomes      = flx.ListProp(  [] , settable=True, doc="List of chromosomes")
    samples          = flx.ListProp(  [] , settable=True, doc="List of samples")

    genome           = flx.StringProp("-", settable=True, doc="Genome name")
    bin_width        = flx.StringProp("-", settable=True, doc="Bin Width name")
    metric           = flx.StringProp("-", settable=True, doc="Metric name")
    chromosome       = flx.StringProp("-", settable=True, doc="Genome name")
    sample           = flx.StringProp("-", settable=True, doc="Sample name")

    genome_data      = flx.DictProp(  {} , settable=True, doc="Genome info")
    chromosome_data  = flx.DictProp(  {} , settable=True, doc="Chromosome info")

    _sel_genomes     = "Select genome"
    _sel_bin_widths  = "Select bin width"
    _sel_metrics     = "Select metric"
    _sel_chromosomes = "Select chromosome"
    _sel_samples     = "Select sample"

    color            = flx.StringProp(color_maps[list(color_maps.keys())[0]][0], settable=True , doc="Color name")
    colors_data      = flx.DictProp(color_maps, settable=False, doc="List of colors")
    _sel_colors      = "Select color"


    def _getElementById(self, eid):
        global window
        el = window.document.getElementById(eid)
        return el

    def _get_and_set(self, eid, f):
        # print("Forms._get_and_set", eid, f)
        el = self._getElementById(eid)
        # print("Forms._get_and_set", eid, f, el)
        if el:
            val = el.value
            # print("Forms._get_and_set", eid, f, el, val)
            if val:
                f(val)

    def _create_dom(self):
        print("Forms._create_dom")
        # Use this method to create a root element for this widget.
        # If you just want a <div> you don't have to implement this.

        main_div = flx.create_element('div', {"id": "main", "class": "selector"}, "") # the default is <div>
        # self.reaction
        return main_div

    def _download_image(self, *ev):
        print("_download_image")
        # var canvas = document.getElementById("mycanvas");
        #     var img    = canvas.toDataURL("image/png");
        #     document.write('<img src="'+img+'"/>');

    def _render_dom(self):
        print("Forms._render_dom")
        # Use this to determine the content. This method may return a
        # string, a list of virtual nodes, or a single virtual node
        # (which must match the type produced in _create_dom()).
        # genomes = self.root.genomes.genomes()

        _getElementById = self._getElementById
        _get_and_set    = self._get_and_set
        _cel            = flx.create_element

        els = []
                    
        def gen_opts(options, option, place_holder_name, place_holder_value, class_str, element_id, callback):
            # print("options", options, "option", option, "place_holder_name", place_holder_name, "place_holder_value", place_holder_value, "class_str", class_str, "element_id", element_id, "callback", callback)
            lopts  = [_cel('option', {"value": place_holder_value}, place_holder_name)]
            lopts += [_cel('option', {"value": s, "selected": "selected" if s == option else None}, s) for s in options]
            lel    =  _cel('select', {"class": class_str, "disabled": False if options else True, "id": element_id, "onchange": gen_getter_and_setter(element_id, callback)}, *lopts)
            return lel

        def gen_getter_and_setter(eid, f):
            return lambda x: _get_and_set(eid, f)

        self.sel_genomes_el     = gen_opts(self.genomes    , self.genome    , self._sel_genomes    , "-", "form_sel", "sel_genomes"    , self.set_genome    )
        self.sel_bin_widths_el  = gen_opts(self.bin_widths , self.bin_width , self._sel_bin_widths , "-", "form_sel", "sel_bin_widths" , self.set_bin_width )
        self.sel_metrics_el     = gen_opts(self.metrics    , self.metric    , self._sel_metrics    , "-", "form_sel", "sel_metrics"    , self.set_metric    )
        self.sel_chromosomes_el = gen_opts(self.chromosomes, self.chromosome, self._sel_chromosomes, "-", "form_sel", "sel_chromosomes", self.set_chromosome)
        self.sel_samples_el     = gen_opts(self.samples    , self.sample    , self._sel_samples    , "-", "form_sel", "sel_samples"    , self.set_sample    )

        els += [
            self.sel_genomes_el,
            self.sel_bin_widths_el,
            self.sel_metrics_el,
            self.sel_chromosomes_el,
            self.sel_samples_el
        ]

        if self.colors_data:
            colors_ph    = _cel('option', {"value": "-"}, self._sel_colors)
            colors_opts  = [colors_ph]

            for color_group, color_names in self.colors_data.items():
                color_names_opts = []
                for color_name in color_names:
                    color_name_opt = _cel('option', {"value": color_name, "selected": "selected" if color_name == self.color else None}, color_name)
                    color_names_opts.append(color_name_opt)
                color_opt = _cel('optgroup', {"label": color_group}, *color_names_opts)
                colors_opts.append(color_opt)
                # colors_opts.extend(color_names_opts)
            self.sel_colors_el = _cel('select', {"id": "sel_colors", "class": "form_sel", "onchange": gen_getter_and_setter("sel_colors", self.set_color)}, *colors_opts)
            els += [self.sel_colors_el]

        _download_image    = self._download_image
        btn_download_image = _cel('button', {"class": "btn btn_download", "onclick": _download_image}, "Download")
        els += [btn_download_image]

        # for k,v in [
        #         ["Genome"    , self.genome    ],
        #         ['Bin Width' , self.bin_width ],
        #         ['Metric'    , self.metric    ],
        #         ['Chromosome', self.chromosome],
        #         ['Sample'    , self.sample    ],
        #         ['Color'     , self.color     ],
        #     ]:
        #     els.append(_cel('span', {'class': 'info_selection_sp'}, _cel('span', {'class': 'info_selection_he'}, k +':'), _cel('span', {'class': 'info_selection_va'}, v )))

        if  self.genome_data:
            for k,v in self.genome_data.items():
                els.append(
                    _cel('span', {'class': 'info_genome_sp info'},
                        _cel('span', {'class': 'info_genome_he'}, " ".join(k.split("_")).title() +':'),
                        _cel('span', {'class': 'info_genome_va'}, str(v) )))

        els += [_cel('br')]

        if self.chromosome_data:
            for k,v in self.chromosome_data.items():
                els.append(
                    _cel('span', {'class': 'info_chromosome_sp info'},
                        _cel('span', {'class': 'info_chromosome_he'}, " ".join(k.split("_")).title() +':'),
                        _cel('span', {'class': 'info_chromosome_va'}, str(v) )))

        els += [_cel('br')]

        return els



    @flx.action
    def reset_genomes(self):
        if self.genomes:
            print("Forms.reset_genomes")
            self.set_genomes([])

    @flx.action
    def reset_bin_widths(self):
        if self.bin_widths:
            print("Forms.reset_bin_widths")
            self.set_bin_widths([])

    @flx.action
    def reset_metrics(self):
        if self.metrics:
            print("Forms.reset_metrics")
            self.set_metrics([])

    @flx.action
    def reset_chromosomes(self):
        if self.chromosomes:
            print("Forms.reset_chromosomes")
            self.set_chromosomes([])

    @flx.action
    def reset_samples(self):
        if self.samples:
            print("Forms.reset_samples")
            self.set_samples([])

    @flx.action
    def reset_canvas(self):
        print("Forms.reset_canvas")



    @flx.action
    def reset_genome_data(self):
        if self.genome_data:
            print("Forms.reset_genome_data")
            self.set_genome_data({})

    @flx.action
    def reset_chromosome_data(self):
        if self.chromosome_data:
            print("Forms.reset_chromosome_data")
            self.set_chromosome_data({})



    @flx.action
    def reset_genome(self):
        if self.genome != "-":
            print("Forms.reset_genome")
            self.set_genome("-")

    @flx.action
    def reset_bin_width(self):
        if self.bin_width != "-":
            print("Forms.reset_bin_width")
            self.set_bin_width("-")

    @flx.action
    def reset_metric(self):
        if self.metric != "-":
            print("Forms.reset_metric")
            self.set_metric("-")

    @flx.action
    def reset_chromosome(self):
        if self.chromosome != "-":
            print("Forms.reset_chromosome")
            self.set_chromosome("-")

    @flx.action
    def reset_sample(self):
        if self.sample != "-":
            print("Forms.reset_sample")
            self.set_sample("-")



    @flx.reaction('genomes')
    def reaction_set_genomes(self, *ev):
        print("Forms.reaction_set_genomes", JSON.stringify(self.genomes))
        self.reset_genome()
        self.reset_bin_widths()
        if self.genomes:
            self.set_genome(self.genomes[0])

    @flx.reaction('bin_widths')
    def reaction_set_bin_widths(self, *ev):
        print("Forms.reaction_set_bin_widths", JSON.stringify(self.bin_widths))
        self.reset_bin_width()
        self.reset_metrics()
        if self.bin_widths:
            self.set_bin_width(self.bin_widths[0])

    @flx.reaction('metrics')
    def reaction_set_metrics(self, *ev):
        print("Forms.reaction_set_metrics", JSON.stringify(self.metrics))
        self.reset_metric()
        self.reset_chromosomes()
        if self.metrics:
            self.set_metric(self.metrics[0])

    @flx.reaction('chromosomes')
    def reaction_set_chromosomes(self, *ev):
        print("Forms.reaction_set_chromosomes", JSON.stringify(self.chromosomes))
        self.reset_chromosome()
        self.reset_chromosome_data()
        self.reset_samples()
        if self.chromosomes:
            self.set_chromosome(self.chromosomes[0])

    @flx.reaction('samples')
    def reaction_set_samples(self, *ev):
        print("Forms.reaction_set_samples")#, JSON.stringify(self.samples))
        self.reset_sample()
        self.reset_canvas()
        if self.samples:
            self.set_sample(self.samples[0])

    @flx.reaction('colors_data')
    def reaction_set_colors_data(self, *ev):
        print("Forms.reaction_colors_data", JSON.stringify(self.colors_data))
        self.reset_canvas()
        if self.colors_data:
            color = self.colors_data[list(self.colors_data.keys())[0]][0]
            print("SETTING COLOR", color)
            self.set_color(color)



    @flx.reaction('genome')
    def reaction_set_genome(self, *ev):
        print("Forms.reaction_set_genome", self.genome)
        self.reset_bin_widths()
        if self.genome != "-":
            # self.reset_bin_width()
            # self.root.update_bin_widths()

            self.root.ui_update_genome_info(self.genome)
            self.root.ui_update_bin_widths(self.genome)

    @flx.reaction('bin_width')
    def reaction_set_bin_width(self, *ev):
        print("Forms.reaction_set_bin_width", self.bin_width)
        # self.reset_metric()
        self.reset_metrics()
        # self.root.update_metrics()
        if self.bin_width != "-":
            self.root.ui_update_metrics(self.genome, self.bin_width)

    @flx.reaction('metric')
    def reaction_set_metric(self, *ev):
        print("Forms.reaction_set_metric", self.metric)
        # self.reset_chromosome()
        self.reset_chromosomes()
        # self.root.update_chromosomes()
        if self.metric != "-":
            self.root.ui_update_chromosomes(self.genome, self.bin_width, self.metric)

    @flx.reaction('chromosome')
    def reaction_set_chromosome(self, *ev):
        print("Forms.reaction_set_chromosome", self.chromosome)
        # self.reset_sample()
        self.reset_samples()
        # self.root.update_samples()
        if self.chromosome != "-":
            self.root.ui_update_chromosome_info(self.genome, self.bin_width, self.metric, self.chromosome)
            self.root.ui_update_samples(self.genome, self.bin_width, self.metric, self.chromosome)

    @flx.reaction('sample')
    def reaction_set_sample(self, *ev):
        print("Forms.reaction_set_sample", self.sample)
        if self.sample != "-" and self.color != "-":
            print("Forms.reaction_set_sample", self.sample, "UPDATING CANVAS")
            self.root.ui_update_canvas(self.genome, self.bin_width, self.metric, self.chromosome, self.sample, self.color)

    @flx.reaction('color')
    def reaction_set_color(self, *ev):
        print("Forms.reaction_set_color", self.color)
        if self.sample != "-" and self.color != "-":
            print("Forms.reaction_set_color", self.sample, "UPDATING CANVAS")
            self.root.ui_update_canvas(self.genome, self.bin_width, self.metric, self.chromosome, self.sample, self.color)



class ChromosomeController(flx.PyComponent):
    def init(self):
        super().init()

        self.chromosome                : reader.Chromosome = None
        self.filename                  : str = "-"

        self.vcf_name                  : str = "-"
        self.bin_width                 : int = -1
        self.chromosome_order          : int = -1
        self.chromosome_name           : str = "-"
        self.metric                    : str = "-"

        self.matrix_size               : int = -1
        self.bin_max                   : int = -1
        self.bin_min                   : int = -1
        self.bin_count                 : int = -1

        self.bin_snps_min              : int = -1
        self.bin_snps_max              : int = -1
        self.bin_snps                  : int = []

        self.chromosome_snps           : int = -1
        self.chromosome_first_position : int = -1
        self.chromosome_last_position  : int = -1

        self.sample_names              : typing.List[str] = []
        self.sample_count              : int = -1
        self.sample_name               : str = "-"

        self.min_val                   : float = 0.0
        self.max_val                   : float = 1.0
        self.num_vals                  : int   = 256

        self.color_name                : str = "Greys"
        self.color_bad                 : str = "pink"
        self.color_over                : str = "pink"
        self.color_under               : str = "pink"

        self.cmap   = None
        self.acolor = None

        self.reset()

    def reset(self):
        print("ChromosomeController.reset")
        self.chromosome                = None

        self.filename                  = "-"
    
        self.vcf_name                  = "-"
        self.bin_width                 = -1
        self.chromosome_order          = -1
        self.chromosome_name           = "-"
        self.metric                    = "-"

        self.matrix_size               = -1
        self.bin_max                   = -1
        self.bin_min                   = -1
        self.bin_count                 = -1

        self.bin_snps_min              = -1
        self.bin_snps_max              = -1
        self.bin_snps                  = []

        self.chromosome_snps           = -1
        self.chromosome_first_position = -1
        self.chromosome_last_position  = -1

        self.sample_names              = []
        self.sample_count              = -1
        self.sample_name               = "-"

        print("ChromosomeController.reset :: done")

    def update(self):
        print("ChromosomeController.update")

        if self.chromosome is None:
            print("ChromosomeController.update :: NO CHROMOSOME")

        else:
            print("ChromosomeController.update :: updating")
            self.file_name                 = self.chromosome.file_name

            self.vcf_name                  = self.chromosome.vcf_name
            self.bin_width                 = self.chromosome.bin_width
            self.chromosome_order          = self.chromosome.chromosome_order
            self.chromosome_name           = self.chromosome.chromosome_name
            self.metric                    = self.chromosome.metric

            self.matrix_size               = self.chromosome.matrix_size
            self.bin_max                   = self.chromosome.bin_max
            self.bin_min                   = self.chromosome.bin_min
            self.bin_count                 = self.chromosome.bin_count

            self.bin_snps_min              = self.chromosome.bin_snps_min
            self.bin_snps_max              = self.chromosome.bin_snps_max
            self.bin_snps                  = self.chromosome.binsnpNp.tolist()

            self.chromosome_snps           = self.chromosome.chromosome_snps
            self.chromosome_first_position = self.chromosome.chromosome_first_position
            self.chromosome_last_position  = self.chromosome.chromosome_last_position

            self.sample_names              = self.chromosome.sample_names
            self.sample_count              = self.chromosome.sample_count
            self.sample_name               = "-"

    def set_chromosome(self, chromosome: reader.Chromosome):
        print("ChromosomeController.set_chromosome")

        self.reset()

        self.chromosome = chromosome

        if self.chromosome is not None:
            self.update()

    def set_min_val(self, min_val: float):
        print("ChromosomeController.set_min_val", min_val)
        self.min_val           = min_val
        self.cmap, self.acolor = self.gen_color()

    def set_max_val(self, max_val: float):
        print("ChromosomeController.set_max_val", max_val)
        self.max_val           = max_val
        self.cmap, self.acolor = self.gen_color()

    def set_min_max_val(self, min_val: float, max_val: float):
        print("ChromosomeController.set_min_max_val", min_val, max_val)
        self.min_val = min_val
        self.max_val = max_val
        self.cmap, self.acolor = self.gen_color()

    def set_sample_name(self, sample_name: str, display=True):
        print("ChromosomeController.set_sample_name", self.sample_name)
        self.sample_name = sample_name
        if display:
            self.display()

    def set_color(self, color_name=None, min_val=None, max_val=None, display=True):
        if  (
                (color_name is not None and self.color_name != color_name) or
                (min_val    is not None and self.min_val    != min_val   ) or
                (max_val    is not None and self.max_val    != max_val   )
        ):
            self.color_name = self.color_name if color_name is None else color_name
            self.min_val    = self.min_val    if min_val    is None else min_val
            self.max_val    = self.max_val    if max_val    is None else max_val

            self.cmap, self.acolor = self.gen_color()
            if display:
                self.display()

    def gen_color(self, color_name=None, min_val=None, max_val=None):
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
        # self.cmap   = cmap
        # self.acolor = acolor
        return cmap, acolor

    def display(self):
        print(f"ChromosomeController.display :: sample_name {self.sample_name} color_name {self.color_name}")

        if self.sample_name == "-":
            print(f"ChromosomeController.display :: sample_name {self.sample_name} color_name {self.color_name} - INVALID SAMPLE NAME")
            return

        if self.color_name == "-":
            print(f"ChromosomeController.display :: sample_name {self.sample_name} color_name {self.color_name} - INVALID COLOR NAME")
            return

        if not (
            self.chromosome      is not None and
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
            bin_cmap, bin_acolor = self.gen_color(min_val=bin_snps_min, max_val=bin_snps_max)
            bin_cmapg  = bin_cmap(self.bin_snps)
            bin_cmapg  = bin_cmapg.tolist()
            bin_acolor = bin_acolor.tolist()
            # print("bin_cmapg", bin_cmapg)

            # print(f"ChromosomeController.display :: matrix = {matrix}")
            print(f"ChromosomeController.display :: matrix max = {matrix_max}")
            print(f"ChromosomeController.display :: matrix min = {matrix_min}")
            self.set_min_val(matrix_min)
            self.set_max_val(matrix_max)

            self.set_color(min_val=matrix_min, max_val=matrix_max)
            img                    = self.cmap(matrix)

            # print(f"ChromosomeController.display :: img = {img}")
            print(f"ChromosomeController.display :: img shape  = {img.shape}")
            print(f"ChromosomeController.display :: img type   = {type(img)}")
            print(f"ChromosomeController.display :: bin_count  = {self.bin_count}")
            # print(f"ChromosomeController.display :: bin_count  = {self.bin_n}")
            imgd      = img.tolist()
            acolord   = self.acolor.tolist()
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

            print(f"ChromosomeController.display :: triggering")

            self.root.graph.set_points({
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
            print(f"ChromosomeController.display :: done")

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



class MainController(flx.PyWidget):
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

        self.genomes    = MainController.genomes_inst
        self.chromosome = ChromosomeController()
        self.forms      = Forms()
        self.graph      = Graph()
        self.verbose    = False
        self.genomes_update()

        self.ui_update_genomes()
        # self.BrowserUI.update_genome_names()

    @flx.action
    def ui_update_genomes(self):
        print("MainController.ui_update_genomes")

        genome_names = self.genomes_names()

        print("MainController.ui_update_genomes", genome_names)
        self.forms.set_genomes(genome_names)

    @flx.action
    def ui_update_genome_info(self, genome_name: str):
        print("MainController.ui_update_genome_info :: genome_name", genome_name)

        genome_info   = self.genomes_genome_info(genome_name)
        genome_info_f = {k:v for k,v in genome_info.items() if isinstance(v, (int, float, str))}

        # print(list(genome_info_f.keys()))
        genome_info_f.update({k:os.path.basename(genome_info_f[k]) for k in ["project_path"]})

        # print("MainController.ui_update_genome_info :: genome_info", genome_info_f)
        self.forms.set_genome_data(genome_info_f)

    @flx.action
    def ui_update_bin_widths(self, genome_name: str):
        print("MainController.ui_update_bin_widths :: genome_name", genome_name)

        bin_widths   = self.genomes.bin_widths(genome_name)
        bin_widths_f = [str(b) for b in bin_widths]

        print("MainController.ui_update_bin_widths :: bin_widths", bin_widths_f)
        self.forms.set_bin_widths(bin_widths_f)

    @flx.action
    def ui_update_metrics(self, genome_name: str, bin_width: str):
        print("MainController.ui_update_metrics :: genome_name", genome_name, "bin_width", bin_width)
        bin_width = int(bin_width)

        metrics   = self.genomes_metrics(genome_name, bin_width)

        print("MainController.ui_update_metrics :: metrics", metrics)
        self.forms.set_metrics(metrics)
    
    @flx.action
    def ui_update_chromosomes(self, genome_name: str, bin_width: str, metric: str):
        print("MainController.ui_update_chromosomes :: genome_name", genome_name, "bin_width", bin_width, "metric", metric)
        bin_width        = int(bin_width)

        chromosome_names = self.genomes_chromosome_names(genome_name, bin_width, metric)
        
        # print("MainController.ui_update_chromosomes :: chromosome_names", chromosome_names)
        self.forms.set_chromosomes(chromosome_names)
    
    @flx.action
    def ui_update_chromosome_info(self, genome_name: str, bin_width: str, metric: str, chromosome: str):
        print("MainController.ui_update_chromosome_info :: genome_name", genome_name, "bin_width", bin_width, "metric", metric, "chromosome", chromosome)
        bin_width         = int(bin_width)

        chromosome_info   = self.chromosome_info(genome_name, bin_width, metric, chromosome)

        chromosome_info_f = {k:v for k,v in chromosome_info.items() if k not in ["sample_names"] and not k.startswith("type_")}
        chromosome_info_f.update({k:os.path.basename(chromosome_info_f[k]) for k in ["vcf_name"]})

        # print(list(chromosome_info_f.keys()))

        # print("MainController.ui_update_chromosome_info :: chromosome_info", chromosome_info_f)
        self.forms.set_chromosome_data(chromosome_info_f)

    @flx.action
    def ui_update_samples(self, genome_name: str, bin_width: str, metric: str, chromosome: str):
        print("MainController.ui_update_samples :: genome_name", genome_name, "bin_width", bin_width, "metric", metric, "chromosome", chromosome)
        bin_width = int(bin_width)

        samples = self.genome_samples(genome_name, bin_width, metric, chromosome)

        # print("MainController.ui_update_samples :: samples", samples)
        self.forms.set_samples(samples)
    
    @flx.action
    def ui_update_canvas(self, genome_name: str, bin_width: str, metric: str, chromosome: str, sample_name: str, color_name: str):
        print("MainController.ui_update_canvas :: genome_name", genome_name, "bin_width", bin_width, "metric", metric, "chromosome", chromosome, "sample_name", sample_name, "color_name", color_name)
        bin_width = int(bin_width)

        self.genomes_load_chromosome(genome_name, bin_width, metric, chromosome)
        self.chromosome.set_color(color_name, display=False)
        self.chromosome.set_sample_name(sample_name, display=False)
        self.chromosome.display()

    ##
    ## Genomes
    ##
    def genomes_reset(self):
        print("MainController.genomes_reset")
        self.genomes = None
        self.genomes.genome.reset()

    def genomes_set(self, genomes: reader.Genomes):
        print(f"MainController.genomes_set")
        self.genomes_reset()
        self.genomes = genomes

    # def genomes_update_genome(self, genome_name: str, bin_width: int, metric: str):
    #     print(f"MainController.genomes_update_genome :: genome_name {genome_name} bin_width {bin_width} metric {metric}")
    #     self.load_genome(genome_name, bin_width, metric)
    #     # self.genome.set_genome(genome_name, genome)

    # def genomes_update_chromosome(self, genome_name: str, bin_width: int, metric: str, chromosome_name: str):
    #     print(f"MainController.genomes_update_chromosome :: genome_name {genome_name} bin_width {bin_width} metric {metric} chromosome_name {chromosome_name}")
    #     self.genomes_load_genome(genome_name, bin_width, metric, chromosome_name)

    def genomes_names(self):
        print(f"MainController.genomes_names")
        return self.genomes.genomes

    def genomes_genome_info(self, genome_name: str):
        print(f"MainController.genomes_genome_info")
        return self.genomes.genome_info(genome_name)

    def genomes_bin_widths(self, genome_name: str) -> typing.List[str]:
        print(f"MainController.genomes_bin_widths")
        return self.genomes.bin_widths(genome_name)

    def genomes_bin_width_info(self, genome_name: str, bin_width: int):
        print(f"MainController.genomes_bin_width_info")
        return self.genomes.bin_width_info(genome_name, bin_width)

    def genomes_metrics(self, genome_name: str, bin_width: int) -> typing.List[str]:
        print(f"MainController.genomes_metrics")
        return self.genomes.metrics(genome_name, bin_width)

    def genomes_metric_info(self, genome_name: str, bin_width: int, metric: str) -> typing.List[str]:
        print(f"MainController.genomes_metric_info")
        return self.genomes.metric_info(genome_name, bin_width, metric)

    def genomes_chromosome_names(self, genome_name: str, bin_width: int, metric: str) -> typing.List[typing.Tuple[int, str]]:
        print(f"MainController.genomes_chromosome_names")
        val = self.genomes.chromosome_names(genome_name, bin_width, metric)
        val = [v[1] for v in val]
        return val

    def genomes_update(self):
        print(f"MainController.genomes_update")
        self.genomes.update(verbose=self.verbose)

    def genomes_load_genome(self, genome_name: str, bin_width: int, metric: str) -> reader.Genome:
        print(f"MainController.genomes_load_genome")
        return self.genomes.load_genome(genome_name, bin_width, metric)

    def genomes_load_chromosome(self, genome_name: str, bin_width: int, metric: str, chromosome_name: str) -> reader.Chromosome:
        print(f"MainController.genomes_load_chromosome")
        chromosome = self.genomes.load_chromosome(genome_name, bin_width, metric, chromosome_name)
        self.chromosome.set_chromosome(chromosome)
        return chromosome

    def genome_samples(self, genome_name: str, bin_width: int, metric: str, chromosome_name: str) -> typing.List[str]:
        print(f"MainController.genome_samples")
        chromosome = self.genomes_load_chromosome(genome_name, bin_width, metric, chromosome_name)
        return chromosome.sample_names

    def chromosome_info(self, genome_name: str, bin_width: int, metric: str, chromosome_name: str):
        print(f"MainController.chromosome_info")
        chromosome = self.genomes_load_chromosome(genome_name, bin_width, metric, chromosome_name)
        return chromosome.metadata()



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


"""
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
            val = self.root.genomes.genome_names()
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
"""
