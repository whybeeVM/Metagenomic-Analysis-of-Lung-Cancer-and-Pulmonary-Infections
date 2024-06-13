#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import argparse
import logging
from collections import defaultdict


LOG = logging.getLogger(__name__)

__version__ = "20211026"
__author__ = ("Liu huifang",)
__email__ = "hfliu@visionmedicals.com"
__all__ = []


def merge_fc(args):
    input_file = dict()
    input_dir = args.dir

    for line in os.listdir(input_dir):
        if line.endswith(".hg38_count.txt"):
            tem = line.strip().split("/")
            tem[-1] = tem[-1].replace(".hg38_count.txt","")
            input_file[tem[-1]] = line

    all_feature = defaultdict(list)
    for k,v in input_file.items():
        order = []
        for sub in open("%s/%s" % (input_dir, v), 'r'):
            if sub.startswith("#"):
                continue
            subtemp = sub.strip().split()
            order.append(subtemp[0])
            if subtemp[0] == "Geneid":
                all_feature[subtemp[0]].append(k)
            else:
                tmp = int(subtemp[6])
                all_feature[subtemp[0]].append(str(tmp))

    f_out = open(args.output, 'w')
    for k in order:
        info_line = '\t'.join(all_feature[k])
        f_out.write('%s\t%s\n' % (k, info_line))
    f_out.close()


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
        merge featurecount files
version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    parser.add_argument('-d',"--dir", required=True,
                        help='the dir of featurecount result')
    parser.add_argument('-o', "--output", default="RNAseq_featurecount.txt",
                        help='output file, default=RNAseq_featurecount.txt')
    args = parser.parse_args()

    merge_fc(args)


if __name__ == "__main__":
    main()
