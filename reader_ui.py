#!/usr/bin/env python3

import os
import sys
import typing

import numpy as np

import flexx
from   flexx import flx, ui

import reader

DEBUG       = True



class browserUI(flx.PyWidget):
    genome_name      = flx.StringProp("", settable=True, doc="current genome")
    bin_width        = flx.IntProp(   -1, settable=True, doc="current bin width")
    metric           = flx.StringProp("", settable=True, doc="current metric")
    chromosome_name  = flx.StringProp("", settable=True, doc="current chromosome")

    _combo_genome_names_text     = "Select genome"
    _combo_bin_widths_text       = "Select bin width"
    _combo_metrics_text          = "Select distance metric"
    _combo_chromosome_names_text = "Select chromosome"

    _is_loaded = False

    def init(self, parent):
        super().init()

        self._is_loaded = False

        self.combo_genome_names     = ui.ComboBox(text=self._combo_genome_names_text    , editable=False, options=[], flex=1, parent=parent, css_class="combo_sel")
        self.combo_bin_widths       = ui.ComboBox(text=self._combo_bin_widths_text      , editable=False, options=[], flex=1, parent=parent, css_class="combo_sel")
        self.combo_metrics          = ui.ComboBox(text=self._combo_metrics_text         , editable=False, options=[], flex=1, parent=parent, css_class="combo_sel")
        self.combo_chromosome_names = ui.ComboBox(text=self._combo_chromosome_names_text, editable=False, options=[], flex=1, parent=parent, css_class="combo_sel")

        self._is_loaded = True

        self.setHidden(self.combo_bin_widths      , True)
        self.setHidden(self.combo_metrics         , True)
        self.setHidden(self.combo_chromosome_names, True)

        # self.root.genomes.update_genome_names()

    def setHidden(self, box, status):
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


    @flx.action
    def set_genome_name(self, genome_name: str):
        if self._is_loaded and genome_name != "":
            print("browserUI.set_genome_name", genome_name)
            self._mutate_genome_name(genome_name)
            self.reset_bin_widths()
            self.update_bin_widths()
            self.setHidden(self.combo_bin_widths, False)

    @flx.action
    def set_bin_width(self, bin_width: int):
        if self._is_loaded and bin_width != -1 and bin_width != "":
            print("browserUI.set_bin_width", bin_width)
            self._mutate_bin_width(bin_width)
            self.reset_metrics()
            self.update_metrics()
            self.setHidden(self.combo_metrics, False)

    @flx.action
    def set_metric(self, metric: str):
        if self._is_loaded and metric != "":
            print("browserUI.set_metric", metric)
            self._mutate_metric(metric)
            self.reset_chromosome_names()
            self.update_chromosome_names()
            # self.root.genomes.update_chromosome_names(self.genome_name, self.bin_width, self.metric)

            if (
                self._is_loaded            and
                self.genome_name     != "" and
                self.bin_width       != -1 and self.bin_width != "" and
                self.metric          != ""
            ):
                print("browserUI.set_metric", self.genome_name, self.bin_width, self.metric)

                self.root.genomes.update_genome(self.genome_name, self.bin_width, self.metric)

    @flx.action
    def set_chromosome_name(self, chromosome_name: str):
        if self._is_loaded and chromosome_name != "":
            print("browserUI.set_chromosome_name", chromosome_name)
            self._mutate_chromosome_name(chromosome_name)

            if (
                self._is_loaded and
                self.genome_name     != "" and
                self.bin_width       != -1 and self.bin_width != "" and
                self.metric          != "" and
                self.chromosome_name != ""
            ):
                print("browserUI.set_chromosome_name", self.genome_name, self.bin_width, self.metric, self.chromosome_name)

                self.root.genomes.update_chromosome(self.genome_name, self.bin_width, self.metric, self.chromosome_name)



    # @flx.action
    # def update_genome(self, genome_name: str, bin_width: int, metric: str, genome: reader.Genome):
    #     assert self.genome_name == genome_name
    #     assert self.bin_width   == bin_width
    #     assert self.metric      == metric
    #     print("browserUI.update_genome", genome_name, bin_width, metric)
    #     # self.root.genomes.update_genome(genome_name, bin_width, metric, genome)
    #     # flx.CheckBox(parent=self.box)

    # @flx.action
    # def update_chromosome(self, genome_name: str, bin_width: int, metric: str, chromosome_name: str, chromosomeL: ChromosomeController):
    #     assert self.genome_name     == genome_name
    #     assert self.bin_width       == bin_width
    #     assert self.metric          == metric
    #     assert self.chromosome_name == chromosome_name
    #     print("browserUI.update_chromosome", genome_name, bin_width, metric, chromosome_name)
        # self.root.genomes.update_genome(genome_name, bin_width, metric, genome)



    @flx.action
    def reset_genome_names(self):
        self.setHidden(self.combo_genome_names, True)
        self.set_genome_names([])
        self.set_genome_names("")
        self.combo_genome_names.set_text(self._combo_genome_names_text)
        self.reset_bin_widths()

    @flx.action
    def reset_bin_widths(self):
        self.setHidden(self.combo_bin_widths, True)
        self.set_bin_widths([])
        self.set_bin_width(-1)
        self.combo_bin_widths.set_text(self._combo_bin_widths_text)
        self.reset_metrics()

    @flx.action
    def reset_metrics(self):
        self.setHidden(self.combo_metrics, True)
        self.set_metrics([])
        self.set_metric("")
        self.combo_metrics.set_text(self._combo_metrics_text)
        self.reset_chromosome_names()

    @flx.action
    def reset_chromosome_names(self):
        self.setHidden(self.combo_chromosome_names, True)
        self.set_chromosome_names([])
        self.set_chromosome_name("")
        self.combo_chromosome_names.set_text(self._combo_chromosome_names_text)



    @flx.action
    def update_genome_names(self):
        print("browserUI.update_genome_names")
        val = self.root.genomes.genomes()
        self.set_genome_names(val)

    @flx.action
    def update_bin_widths(self):
        print("browserUI.update_bin_widths")
        val = self.root.genomes.bin_widths(self.genome_name)
        self.set_bin_widths(val)

    @flx.action
    def update_metrics(self):
        print("browserUI.update_metrics")
        val = self.root.genomes.metrics(self.genome_name, self.bin_width)
        self.set_metrics(val)

    @flx.action
    def update_chromosome_names(self):
        print("browserUI.update_chromosome_names")
        val = self.root.genomes.chromosome_names(self.genome_name, self.bin_width, self.metric)
        self.set_chromosome_names(val)



    @flx.action
    def set_genome_names(self, val):
        print("browserUI.set_genome_names")
        self.combo_genome_names.set_options(val)

    @flx.action
    def set_bin_widths(self, val):
        print("browserUI.set_bin_widths")
        self.combo_bin_widths.set_options(val)

    @flx.action
    def set_metrics(self, val):
        print("browserUI.set_metrics")
        self.combo_metrics.set_options(val)

    @flx.action
    def set_chromosome_names(self, val):
        print("browserUI.set_chromosome_names")
        self.combo_chromosome_names.set_options(val)



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

    @flx.reaction("combo_chromosome_names.selected_key")
    def reaction_combo_chromosome_names(self, ev):
        self.set_chromosome_name(self.combo_chromosome_names.selected_key)




    @flx.action
    def update_genome(self):
        print("browserUI.update_genome")
        self.setHidden(self.combo_chromosome_names, False)

    @flx.action
    def update_chromosome(self):
        print("browserUI.update_chromosome")



class ChromosomeController(flx.PyWidget):
    def init(self, parent):
        super().init()

        self.chromosome : reader.Chromosome = None

        flx.Label(text="chromosome", flex=27, parent=parent)

    def set_chromosome(self, chromosome: reader.Chromosome):
        print("ChromosomeController.set_chromosome")
        self.chromosome = chromosome
        self.root.browserUI.update_chromosome()

    @property
    def vcf_name(self) -> str:
        return self.chromosome.vcf_name

    @property
    def bin_width(self) -> int:
        return self.chromosome.bin_width

    @property
    def chromosome_order(self) -> int:
        return self.chromosome.chromosome_order

    @property
    def chromosome_name(self) -> str:
        return self.chromosome.chromosome_name

    @property
    def metric(self) -> str:
        return self.chromosome.metric



    @property
    def matrix_size(self) -> int:
        return self.chromosome.matrix_size

    @property
    def bin_max(self) -> int:
        return self.chromosome.bin_max

    @property
    def bin_min(self) -> int:
        return self.chromosome.bin_min

    @property
    def bin_count(self) -> int:
        return self.chromosome.bin_count



    @property
    def bin_snps_min(self) -> int:
        return self.chromosome.bin_snps_min

    @property
    def bin_snps_max(self) -> int:
        return self.chromosome.bin_snps_max

    @property
    def chromosome_snps(self) -> int:
        return self.chromosome.chromosome_snps

    @property
    def chromosome_first_position(self) -> int:
        return self.chromosome.chromosome_first_position

    @property
    def chromosome_last_position(self) -> int:
        return self.chromosome.chromosome_last_position



    @property
    def sample_names(self) -> reader.SampleNamesType:
        return self.chromosome.sample_names
            
    @property
    def sample_count(self) -> int:
        return self.chromosome.sample_count


    @property
    def filename(self) -> str:
        return self.chromosome.filename

    @property
    def matrix(self) -> np.ndarray:
        return self.chromosome.matrix

    @property
    def matrix_dtype(self) -> np.dtype:
        return self.chromosome.matrix_dtype



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


    # self.matrixNp                     : np.ndarray = None
    # self.matrixNp                     = np.zeros((self.bin_max, self.matrix_size ), self.type_matrix_counter  )

    # self.pairwiNp                     : np.ndarray = None
    # self.pairwiNp                     = np.zeros((self.bin_max, self.sample_count), self.type_pairwise_counter)

    # self.binsnpNp                     : np.ndarray = None
    # self.binsnpNp                     = np.zeros( self.bin_max,                     self.type_pairwise_counter)

    # self.alignmentNp                  : np.ndarray = None
    # self.alignmentNp                  = np.zeros((self.bin_max, self.sample_count), np.unicode_)

    # self.positionNp                   : np.ndarray = None
    # self.positionNp                   = np.zeros((self.bin_max, position_size    ), self.type_positions)



class GenomeController(flx.PyWidget):
    def init(self, parent):
        super().init()

        self.genome     : reader.Genome        = None

        # with flx.HFix(flex=19) as self.chromosomeBox:
        #     self.chromosome = ChromosomeController()
            # ui.Label(text=lambda:"vcf_name {}".format(self.vcf_name))
        # ui.Label(text=lambda:"vcf_name {}".format(self.vcf_name))
        flx.Label(text="genome", flex=1, parent=parent)

    def set_genome(self, genome: reader.Genome):
        print("GenomeController.set_genome")
        self.genome = genome
        self.root.browserUI.update_genome()

    @property
    def vcf_name(self) -> str:
        if self.genome is None:
            return None
        return self.genome.vcf_name

    @property
    def bin_width(self) -> int:
        if self.genome is None:
            return None
        return self.genome.bin_width

    @property
    def metric(self) -> str:
        if self.genome is None:
            return None
        return self.genome.metric

    @property
    def sample_names(self) -> reader.SampleNamesType:
        if self.genome is None:
            return None
        return self.genome.sample_names

    @property
    def sample_count(self) -> int:
        if self.genome is None:
            return None
        return self.genome.sample_count

    @property
    def chromosome_names(self) -> reader.ChromosomeNamesType:
        if self.genome is None:
            return None
        return self.genome.chromosome_names

    @property
    def chromosome_count(self) -> int:
        if self.genome is None:
            return None
        return self.genome.chromosome_count

    @property
    def genome_bins(self) -> int:
        if self.genome is None:
            return None
        return self.genome.genome_bins

    @property
    def genome_snps(self) -> int:
        if self.genome is None:
            return None
        return self.genome.genome_snps

    @property
    def filename(self) -> str:
        if self.genome is None:
            return None
        return self.genome.filename

    # @flx.action
    # def update_genome(self, genome_name: str, bin_width: int, metric: str, genome: reader.Genome):
    #     assert genome_name != "" and genome_name is not None
    #     assert bin_width   != "" and bin_width   is not None and bin_width   != -1
    #     assert metric      != "" and metric      is not None
    #     print("GenomeController.update_genome", genome_name, bin_width, metric)
    #     # self.genome = genome
    #     self.set_genome(genome)

    # def get_chromosome(self, chromosome_name: str):
    #     if self.genome is None:
    #         return None
    #     return ChromosomeController(self.genome.get_chromosome(chromosome_name))



class GenomesController(flx.PyComponent):
    def init(self):
        super().init()

        self._genomes : reader.Genomes   = None

    # def genome(self) -> GenomeController:
    #     return GenomeController(self._genomes.genome)

    # def chromosomes(self) -> reader.ChromosomeNamesType:
    #     return self._genomes.chromosomes

    # def chromosome(self) -> ChromosomeController:
    #     return ChromosomeController(self._genomes.chromosome)

    def set_genomes(self, genomes: reader.Genomes):
        self._genomes = genomes

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



    # @flx.action
    def update_genome(self, genome_name: str, bin_width: int, metric: str):
        print(f"MainController.update_genome genome_name {genome_name} bin_width {bin_width} metric {metric}")
        genome = self.load_genome(genome_name, bin_width, metric)
        self.root.genome.set_genome(genome)

    # @flx.action
    def update_chromosome(self, genome_name: str, bin_width: int, metric: str, chromosome_name: str):
        print(f"MainController.update_chromosome genome_name {genome_name} bin_width {bin_width} metric {metric} chromosome_name {chromosome_name}")
        chromosome = self.load_chromosome(genome_name, bin_width, metric, chromosome_name)
        self.root.chromosome.set_chromosome(chromosome)



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
                self.browserUI  = browserUI(self.header)
            with flx.HFix(flex=1) as self.genomeBox:
                self.genome     = GenomeController(self.genomeBox)
            with flx.HFix(flex=27) as self.chromosomeBox:
                self.chromosome = ChromosomeController(self.chromosomeBox)
            #     self.chromosomeController   = self.root.chromosome

        self.genomes    = GenomesController()

        self.genomes.set_genomes(MainController.genomes_inst)
        self.genomes.update(verbose=True)
        self.browserUI.update_genome_names()



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