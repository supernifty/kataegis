#!/usr/bin/env python
'''
  annotate vcf with inter mutation distance
'''

import argparse
import collections
import logging
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import cyvcf2

DPI=150

COLORS = {
  'C>A': '#00BCF2',
  'C>G': '#000000',
  'C>T': '#E82818',
  'T>A': '#CAC8C9',
  'T>C': '#9BD357',
  'T>G': '#EEC6C4'
}

COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def get_kataegis(v):
  try:
    return v.INFO['kataegis'] == 'YES'
  except KeyError:
    return False

def get_imd(v):
  try:
    return int(v.INFO['imd'])
  except KeyError:
    return None

def color(v):
  if v in COLORS:
    return v, COLORS[v]
  elif len(v) == 3 and v[0] in COMP and v[2] in COMP:
    comp = '{}>{}'.format(COMP[v[0]], COMP[v[2]])
    return comp, COLORS[comp]
  else:
    logging.debug('unrecognized mutation %s', v)
    return 'Other', '#E3FF00'

def main(count, mean, plot_fn, just_kataegis):
  logging.info('reading vcf from stdin...')
  vcf_in = cyvcf2.VCF('-')

  vcf_in.add_info_to_header({'ID': 'imd', 'Description': 'Inter mutation distance', 'Type':'Character', 'Number': '1'})
  vcf_in.add_info_to_header({'ID': 'kataegis', 'Description': 'Inter mutation distance', 'Type':'Character', 'Number': '1'})
  sys.stdout.write(vcf_in.raw_header)

  q = [] # queue of points under consideration
  r = collections.defaultdict(list) # regions of kataegis
  p = collections.defaultdict(list) # points to plot
  m = collections.defaultdict(list) # points to plot
  c = [] # chroms
  kataegis_start = None
  for i, variant in enumerate(vcf_in):
    if variant.CHROM not in c:
      c.append(variant.CHROM)
    if len(q) > 0 and q[-1].CHROM != variant.CHROM:
      if kataegis_start is not None:
        # in kataegis
        r[q[-1].CHROM].append((kataegis_start, kataegis_end))
      for v in q:
        if not just_kataegis or get_kataegis(v):
          sys.stdout.write(str(v))
        if get_imd(v) is not None:
          p[v.CHROM].append((v.POS, '{}>{}'.format(v.REF, v.ALT[0]), get_imd(v), get_kataegis(v)))
          m[v.CHROM].append(None)
      q = []
      kataegis_start = None
      logging.info('processing chromosome %s...', variant.CHROM)

    q.append(variant)
    if len(q) > 1:
      variant.INFO["imd"] = q[-1].POS - q[-2].POS
    if len(q) == count:
      # get all imds in the list
      imds = [get_imd(v) for v in q if get_imd(v) is not None]
      a = sum(imds) / len(imds)
      if a < mean:
        # give them all kataegis
        logging.info('kataegis found with mean %i at %s:%i', a, v.CHROM, v.POS)
        if kataegis_start is None:
          kataegis_start = q[0].POS
        kataegis_end = q[-1].POS
        for v in q:
          v.INFO['kataegis'] = 'YES'
      else: # no kataegis
        #logging.debug('no kataegis found with mean %i at %s:%i', a, variant.CHROM, variant.POS)
        if kataegis_start is not None:
          r[variant.CHROM].append((kataegis_start, kataegis_end))
      # pop item off the queue and write
      v = q.pop(0)
      if not just_kataegis or get_kataegis(v):
        sys.stdout.write(str(v))
      if get_imd(v) is not None:
        p[v.CHROM].append((v.POS, '{}>{}'.format(v.REF, v.ALT[0]), get_imd(v), get_kataegis(v)))
        m[v.CHROM].append(a)

  # write out the queue
  if kataegis_start is not None:
    # in kataegis
    r[q[-1].CHROM].append((kataegis_start, kataegis_end))

  # write out the queue
  for v in q:
    if not just_kataegis or get_kataegis(v):
      sys.stdout.write(str(v))
    if get_imd(v) is not None:
      p[v.CHROM].append((v.POS, '{}>{}'.format(v.REF, v.ALT[0]), get_imd(v), get_kataegis(v)))
      m[v.CHROM].append(None)
 
  # plot
  if plot_fn is not None:
    logging.info('plotting %i points and %i regions of kataegis...', sum([len(p[c]) for c in p]), len(r))
    # points
    xs = []
    ms = []
    ys = []
    cs = []
    ks = []
    chroms = []
    patches = []
    colors = set()
    n = 0
    last_n = 0
    kstart = None
    for ch in c:
      if ch in p:
        for x, mv in zip(p[ch], m[ch]):
          xs.append(n)
          ys.append(max(1, x[2]))
          ms.append(mv)
          if x[3]:
            if kstart is None:
              kstart = n
          else:
            if kstart is not None:
              ks.append((kstart, n))
              kstart = None
          name, col = color(x[1])
          if name not in colors:
            colors.add(name)
            patches.append(mpatches.Patch(color=col, label=name))
          cs.append(col) # todo legend
          n += 1
        if kstart is not None:
          ks.append((kstart, n))
          kstart = None
      chroms.append((ch, last_n, n))
      last_n = n

    logging.info('y range is %i to %i for %i points', min(ys), max(ys), n)
    logging.info('mean range is %i to %i for %i points', min([m for m in ms if m is not None]), max([m for m in ms if m is not None]), n)
    
    #matplotlib.style.use('seaborn')
    fig = plt.figure(figsize=(24, 8))
    ax = fig.add_subplot(111)
    #ax.set_ylim(bottom=1)
    ax.set_yscale('log', nonposy='clip')
    ax.grid(axis='y')
    ax.scatter(xs, ys, c=cs, s=2, alpha=0.5)    
    ax.plot(xs, ms, color='#003090', linewidth=2, alpha=0.5)
    plt.legend(handles=patches)

    # regions
    for k in ks:
      logging.info('kataegis for mutations %i to %i', k[0], k[1])
      ax.axvspan(k[0], k[1], 0, 1, color='#ff0000', alpha=0.3)
  
    # chromosomes
    xticks = []
    xlabels = []
    for c in chroms:
      ax.axvline(c[2], color='#c0c0c0')
      xticks.append(c[1] + (c[2] - c[1])/2)
      xlabels.append(c[0])
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
  
    plt.tight_layout()
    plt.savefig(plot_fn, transparent=False, dpi=DPI)

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--count', default=6, type=int, help='number of mutations to stay below mean')
  parser.add_argument('--mean', default=1000, type=int, help='stay below this mean to identify kataegis')
  parser.add_argument('--plot', required=False, help='plot filename')
  parser.add_argument('--just_kataegis', action='store_true', help='only kataegis is written')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.count, args.mean, args.plot, args.just_kataegis)
