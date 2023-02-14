#install pycirclize
import matplotlib.pyplot as plt
from pycirclize import Circos
from pycirclize.parser import Genbank
import numpy as np
from matplotlib.patches import Patch

# Load Genbank file
gbk_file = ("escherichia_coli.gbk.gz")
gbk = Genbank(gbk_file)

circos = Circos(sectors={gbk.name: gbk.range_size})
circos.text("Escherichia coli\n(NC_000913)", size=12, r=20)
sector = circos.get_sector(gbk.name)

# Plot outer track with xticks
major_ticks_interval = 500000
minor_ticks_interval = 100000
outer_track = sector.add_track((98, 100))
outer_track.axis(fc="lightgrey")
outer_track.xticks_by_interval(
    major_ticks_interval, label_formatter=lambda v: f"{v/ 10 ** 6:.1f} Mb"
)
outer_track.xticks_by_interval(minor_ticks_interval, tick_length=1, show_label=False)

# Plot Forward CDS, Reverse CDS, rRNA, tRNA
f_cds_track = sector.add_track((90, 97), r_pad_ratio=0.1)
f_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=1), fc="red")

r_cds_track = sector.add_track((83, 90), r_pad_ratio=0.1)
r_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=-1), fc="blue")

rrna_track = sector.add_track((76, 83), r_pad_ratio=0.1)
rrna_track.genomic_features(gbk.extract_features("rRNA"), fc="green")

trna_track = sector.add_track((69, 76), r_pad_ratio=0.1)
trna_track.genomic_features(gbk.extract_features("tRNA"), color="magenta", lw=0.1)

# Plot GC content
gc_content_track = sector.add_track((50, 65))

pos_list, gc_contents = gbk.calc_gc_content()
gc_contents = gc_contents - gbk.calc_genome_gc_content()
positive_gc_contents = np.where(gc_contents > 0, gc_contents, 0)
negative_gc_contents = np.where(gc_contents < 0, gc_contents, 0)
abs_max_gc_content = np.max(np.abs(gc_contents))
vmin, vmax = -abs_max_gc_content, abs_max_gc_content
gc_content_track.fill_between(
    pos_list, positive_gc_contents, 0, vmin=vmin, vmax=vmax, color="black"
)
gc_content_track.fill_between(
    pos_list, negative_gc_contents, 0, vmin=vmin, vmax=vmax, color="grey"
)

# Plot GC skew
gc_skew_track = sector.add_track((35, 50))

pos_list, gc_skews = gbk.calc_gc_skew()
positive_gc_skews = np.where(gc_skews > 0, gc_skews, 0)
negative_gc_skews = np.where(gc_skews < 0, gc_skews, 0)
abs_max_gc_skew = np.max(np.abs(gc_skews))
vmin, vmax = -abs_max_gc_skew, abs_max_gc_skew
gc_skew_track.fill_between(
    pos_list, positive_gc_skews, 0, vmin=vmin, vmax=vmax, color="olive"
)
gc_skew_track.fill_between(
    pos_list, negative_gc_skews, 0, vmin=vmin, vmax=vmax, color="purple"
)

fig = circos.plotfig()

# Add legend
handles = [
    Patch(color="red", label="Forward CDS"),
    Patch(color="blue", label="Reverse CDS"),
    Patch(color="green", label="rRNA"),
    Patch(color="magenta", label="tRNA"),
    Patch(color="black", label="Positive GC Content"),
    Patch(color="grey", label="Negative GC Content"),
    Patch(color="olive", label="Positive GC Skew"),
    Patch(color="purple", label="Negative GC Skew"),
]
_ = fig.legend(handles=handles, bbox_to_anchor=(0.5, 0.475), loc="center", fontsize=8)

plt.show()

