#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os

mycolors = [
    None,
    colors.Color(255/255, 68/255, 6/255, alpha=0.5),
    colors.Color(20/255, 166/255, 6/255, alpha=0.5)
]

refrecord = SeqIO.read(os.path.expanduser("~/GitHub/riboSeed/manuscript_results/simulated_genome/mauve/reference.gb"), "genbank")

defererecord = list(SeqIO.parse(os.path.expanduser("~/GitHub/riboSeed/manuscript_results/simulated_genome/de_fere/alignment2/coli_de_fere_novo.fa.fas"), "fasta"))
df_bb = os.path.expanduser("~/GitHub/riboSeed/manuscript_results/simulated_genome/de_fere/alignment2/alignment2.backbone")

denovorecord = list(SeqIO.parse(os.path.expanduser("~/GitHub/riboSeed/manuscript_results/simulated_genome/de_novo/alignment2/coli_de_novo.fa.fas"), "fasta"))

dn_bb = os.path.expanduser("~/GitHub/riboSeed/manuscript_results/simulated_genome/de_novo/alignment2/alignment2.backbone")

bbs = [[], df_bb, dn_bb]

# #Create the feature set and its feature objects,
# gd_feature_set = GenomeDiagram.FeatureSet()
# name = "Proux Fig 6"
# gd_diagram = GenomeDiagram.Diagram(name)

# max_len = 0
# for j, record in enumerate([refrecord, defererecord, denovorecord]):
#     if not isinstance(record, list):
#         max_len = max(max_len, len(record))
#         gd_track_for_features = gd_diagram.new_track(
#             1,
#             name=record.name,
#             greytrack=True,
#             height=.3,
#             start=0, end=len(record))
#         gd_feature_set = gd_track_for_features.new_set()

#         for i, feature in enumerate(record.features):
#             if feature.type != "rRNA":
#                 #Exclude this feature
#                 continue
#             # print (feature)
#             gd_feature_set.add_feature(
#                 feature,
#                 sigil="ARROW",
#                 arrowshaft_height=1,
#                 color=colors.red,
#                 label=True,
#                 name=feature.qualifiers["note"][0],
#                 label_position="start",
#                 label_size=6, label_angle=45)
#     else:
#         this_max = sum([len(x) for x in record])
#         gd_track_for_features = gd_diagram.new_track(
#             j + 1,
#             name="track %s" % str(j + 1),
#             greytrack=True,
#             height=0.2,
#             start=0, end=this_max)
#         gd_feature_set = gd_track_for_features.new_set()
#         pos = 0
#         for contig in record:
#             feature = SeqFeature(FeatureLocation(pos, pos + len(contig)))
#             pos = pos + len(contig)
#             gd_feature_set.add_feature(
#                 feature,
#                 # sigil="OCTO",
#                 # color=gene_colors[i],
#                 color=colors.grey,
#                 border=colors.black,
#                 label=True,
#                 name="lorp",
#                 label_position="start",
#                 label_size=6,
#                 label_angle=45)

# for i, rec in reversed(list(enumerate(comps_list))):
#     # if i != 1:
#     #     continue
#     for line in rec:
#         print(line)
#         track_X = gd_diagram.tracks[1]
#         track_Y = gd_diagram.tracks[i + 1]
#         A_start, A_end, B_start, B_end = line
#         color = mycolors[i]
#         if (B_start == 0 and B_end == 0) or \
#            (A_start == 0 and A_end == 0):
#             continue
#         link_xy = CrossLink((track_X, A_start, A_end),
#                             (track_Y, B_start, B_end),
#                             color,  # fill
#                             border=color)  # border
#         gd_diagram.cross_track_links.append(link_xy)

# print(gd_diagram.__dict__)
# #Create a track, and a diagram
# gd_diagram.draw(format="linear", orientation="landscape", pagesize='A4',
#                 fragments=1, start=0, end=110000)
# gd_diagram.write(os.path.expanduser("~/GitHub/riboSeed/plasmid_linear.pdf"),
#                  "PDF")
from matplotlib.patches import FancyBboxPatch

def draw_bbox(ax, bb):
    # boxstyle=square with pad=0, i.e. bbox itself.
    p_bbox = FancyBboxPatch((bb.xmin, bb.ymin),
                            abs(bb.width), abs(bb.height),
                            boxstyle="square,pad=0.",
                            ec="k", fc="none", zorder=10.,
                            )
    ax.add_patch(p_bbox)


def test1(ax):

    # a fancy box with round corners. pad=0.1
    p_fancy = FancyBboxPatch((bb.xmin, bb.ymin),
                             abs(bb.width), abs(bb.height),
                             boxstyle="round,pad=0.1",
                             fc=(1., .8, 1.),
                             ec=(1., 0.5, 1.))

    ax.add_patch(p_fancy)

    ax.text(0.1, 0.8,
            r' boxstyle="round,pad=0.1"',
            size=10, transform=ax.transAxes)

    # draws control points for the fancy box.
    #l = p_fancy.get_path().vertices
    #ax.plot(l[:,0], l[:,1], ".")

    # draw the original bbox in black
    draw_bbox(ax, bb)


mycolors = {
    "pinkish": mpl.colors.ColorConverter().to_rgba(
        "#ff4c05", alpha=0.2),
    "redish": mpl.colors.ColorConverter().to_rgba(
        "#ff4c05", alpha=0.8),
    "yellish": mpl.colors.ColorConverter().to_rgba(
        "#FFFB07", alpha=0.2),
    "greenish": mpl.colors.ColorConverter().to_rgba(
        "#04FF08", alpha=0.2),
    "bluish": mpl.colors.ColorConverter().to_rgba(
        "#06B9FF", alpha=0.2),
    "greyish": mpl.colors.ColorConverter().to_rgba(
        "#6505FF", alpha=0.2),
    "clear": mpl.colors.ColorConverter().to_rgba(
        "#FF012F", alpha=0.0),
}


def parseBbcols(filelist):
    """ Given a list of .bbcols files, write out as nested list
    """
    comps_list = []
    for i, f in enumerate(filelist):
        if i == 0:
            continue
        with open(f, "r") as infile:
            temp = [x.strip().split("\t") for x in infile.readlines()]
            temp2 = []
            for sublist in temp[1:len(temp)]:
                temp2.append([int(x) for x in sublist])
                # temp = [int(x) for x in [y for y in temp[1:len(temp)]]]
        comps_list.append(temp2)  # get rid of header
    return (comps_list)




def plot_mauve_compare(refgb,
                       assembly_list,
                       bbcols_list,
                       bufferlen=10000,
                       aspect=.6,
                       names=["Position", "Entropy"],
                       title="Shannon Entropy by Position",
                       output_prefix="entropy_plot.png"):
    assert len(assembly_list) == len(bbcols_list), \
        "must have same amount of assemblies as bbcols"
    with open(refgb, "r") as rg:
        ref_recs = list(SeqIO.parse(rg, "genbank"))

    bbcols = parseBbcols(bbcols_list)
    npanels = len(assembly_list) + 1
    ref_combined_len = sum([len(x) for x in ref_recs]) + bufferlen
    print(len(ref_recs[0]))
    print(ref_combined_len)
    fig, axX = plt.subplots(npanels, 1, sharex=True,
                            gridspec_kw={
                                'height_ratios': [1 for x in range(npanels)]})
    axX[0].set_title(title, y=1.08)
    relheight = ref_combined_len * aspect
    relheighteach = relheight / npanels
    xmin, xmax = 0, ref_combined_len
    ymin, ymax = relheighteach * .5, - relheighteach * .5
    for ax in axX:
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
    # yjust = -.1
    # add annotations
    last_chrom_end = 0
    for record in ref_recs:
        # coding sequence
        coding_box = FancyBboxPatch(
            (last_chrom_end, - relheighteach * .05),
            len(record), relheighteach * .1,
            boxstyle="round,pad=0,rounding_size=" + str(relheighteach / 100),
            mutation_aspect=.5,
            # mutation_scale=.5,
            fc=mycolors['greyish'],
            ec=mycolors['clear']
        )
        buffer_box = FancyBboxPatch(
            (last_chrom_end + len(record), - relheighteach * .05),
            last_chrom_end + len(record) + bufferlen, relheighteach * .1,
            boxstyle="round,pad=0,rounding_size=0",
            mutation_aspect=.5,
            # mutation_scale=.5,
            fc=mycolors['clear'],
            ec=mycolors['clear']
        )
        last_chrom_end = last_chrom_end + len(record) + bufferlen
        axX[0].add_patch(coding_box)
        axX[0].add_patch(buffer_box)
        for i, feature in enumerate(record.features):
            if feature.type != "rRNA" and i == 0:
                #Exclude this feature
                continue
            feat_len = \
                feature.location.end.position - feature.location.start.position
            anno_box = FancyBboxPatch(
                (feature.location.start.position, - relheighteach * .1),
                feat_len, relheighteach * .2,
                boxstyle="round,pad=0,rounding_size=" + str(feat_len / 3),
                mutation_aspect=.5,
                # mutation_scale=.5,
                fc=mycolors['redish'],
                ec=mycolors['redish']
            )

            axX[0].add_patch(anno_box)

    for i in range(npanels):
    # for each assembly
        if i == 0:
            continue
        with open(assembly_list[i - 1], "r") as infile:
            contigs = list(SeqIO.parse(infile, "fasta"))
        last_contig_end = 0
        for record in contigs:

            coding_box = FancyBboxPatch(
                (last_contig_end, - relheighteach * .05),
                len(record), relheighteach * .1,
                boxstyle="round,pad=0,rounding_size=" + str(relheighteach / 100),
                mutation_aspect=.5,
                # mutation_scale=.5,
                fc=mycolors['greyish'],
                ec=mycolors['clear']
            )
            buffer_box = FancyBboxPatch(
                (last_contig_end + len(record), - relheighteach * .05),
                last_contig_end + len(record) + bufferlen, relheighteach * .1,
                boxstyle="round,pad=0,rounding_size=0",
                mutation_aspect=.5,
                # mutation_scale=.5,
                fc=mycolors['clear'],
                ec=mycolors['clear']
            )
            last_contig_end = last_contig_end + len(record) + bufferlen
            axX[i].add_patch(coding_box)
            axX[i].add_patch(buffer_box)

    for i, bblist in enumerate(bbcols):
        for block in bblist:

            verts = [
                (0., 0.),  # left, bottom
                (0., 1.),  # left, top
                (1., 1.),  # right, top
                (1., 0.),  # right, bottom
                (0., 0.),  # ignored
            ]

            codes = [mpl.path.Path.MOVETO,
                     mpl.path.Path.LINETO,
                     mpl.path.Path.LINETO,
                     mpl.path.Path.LINETO,
                     mpl.path.Path.CLOSEPOLY]

            path = mpl.path.Path(verts, codes)

            patch = patches.PathPatch(path, facecolor='orange', lw=2)
            ax[i].add_patch(patch)

    # for index, anno in enumerate(anno_list):
    #     rect1 = patches.Rectangle(
    #         (anno[1][0],  # starting x
    #          ymin),  # starting y
    #         anno[1][1] - anno[1][0],  # rel x end
    #         ymax - ymin,  # rel y end
    #         facecolor=mpl.colors.ColorConverter().to_rgba(
    #             colors[index], alpha=0.2),
    #         edgecolor=mpl.colors.ColorConverter().to_rgba(
    #             colors[index], alpha=0.2))
    #     rect2 = patches.Rectangle(
    #         (anno[1][0],  # starting x
    #          1),  # starting y
    #         anno[1][1] - anno[1][0],  # rel x end
    #         cov_max_depth,  # dont -1 beacuse start at 1
    #         facecolor=mpl.colors.ColorConverter().to_rgba(
    #             colors[index], alpha=0.2),
    #         edgecolor=mpl.colors.ColorConverter().to_rgba(
    #             colors[index], alpha=0.2))
    #     ax1.add_patch(rect1)
    #     ax2.add_patch(rect2)
    #     ax1.text((anno[1][0] + anno[1][1]) / 2,    # x location
    #              ymax - 0.48 - yjust,                      # y location
    #              anno[0][0:20],                          # text first 20 char
    #              ha='center', color='red', weight='bold', fontsize=10)
    #     yjust = yjust * - 1
    # ax1.scatter(x=df["Position"], y=df["Entropy"],
    #             marker='o', color='black', s=2)
    # # add smoothing for kicks
    # df["fit"] = savitzky_golay(df["Entropy"].values, 351, 3)  # window size 51, polynomial order 3
    # ax1.scatter(x=df["Position"], y=df["fit"], color='red', s=1)
    # #
    # ax1.set_ylabel('Shannon Entropy')
    # ax1.get_yaxis().set_label_coords(-.05, 0.5)
    # ax2.set_xlim([xmin, xmax])
    # ax2.invert_yaxis()
    # ax2.set_ylabel('Consensus Coverage')
    # ax2.set_xlabel('Position (bp)')
    # ax2.get_yaxis().set_label_coords(-.05, 0.5)
    # # ax2.set_ylim([1, cov_max_depth + 1]) #, 1])
    # ax2.bar(df_con.index, df_con["depth"],
    #         width=1, color='darkgrey', linewidth=0, edgecolor='darkgrey')
    # # ax2.step(df_con.index, df_con["depth"],
    # #          where='mid', color='darkgrey')
    # for ax in [ax1, ax2]:
    #     ax.spines['right'].set_visible(False)
    # # Only show ticks on the left and bottom spines
    # ax1.spines['top'].set_visible(False)
    # ax2.spines['bottom'].set_visible(False)
    # ax.yaxis.set_ticks_position('left')
    # ax2.xaxis.set_ticks_position('bottom')
    # ax1.xaxis.set_ticks_position('top')
    # ax1.tick_params(axis='y', colors='dimgrey')
    # ax2.tick_params(axis='y', colors='dimgrey')
    # ax1.tick_params(axis='x', colors='dimgrey')
    # ax2.tick_params(axis='x', colors='dimgrey')
    # ax1.yaxis.label.set_color('black')
    # ax2.yaxis.label.set_color('black')
    # ax1.xaxis.label.set_color('black')
    # ax2.xaxis.label.set_color('black')
    plt.tight_layout()
    fig.subplots_adjust(hspace=0)
    fig.set_size_inches(12, 7.5)
    fig.savefig(str(output_prefix + '.png'), dpi=(200))
    fig.savefig(str(output_prefix + '.pdf'), dpi=(200))
    return 0

if __name__ == "__main__":
    refgbpath = os.path.expanduser("~/GitHub/riboSeed/manuscript_results/simulated_genome/mauve/reference.gb")

    deferepath = os.path.expanduser("~/GitHub/riboSeed/manuscript_results/simulated_genome/de_fere/alignment2/coli_de_fere_novo.fa.fas")

    df_bb = os.path.expanduser("~/GitHub/riboSeed/manuscript_results/simulated_genome/de_fere/alignment2/alignment2.backbone")

    denovopath = os.path.expanduser("~/GitHub/riboSeed/manuscript_results/simulated_genome/de_novo/alignment2/coli_de_novo.fa.fas")

    dn_bb = os.path.expanduser("~/GitHub/riboSeed/manuscript_results/simulated_genome/de_novo/alignment2/alignment2.backbone")

    plot_mauve_compare(refgb=refgbpath,
                       assembly_list=[deferepath, denovopath],
                       bbcols_list=[df_bb, dn_bb],
                       bufferlen=1000,
                       names=["Position", "Entropy"],
                       title="Shannon Entropy by Position",
                       output_prefix="./test_plot.png")
