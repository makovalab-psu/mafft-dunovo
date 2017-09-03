#!/usr/bin/env python
import os
import sys
import errno
import ctypes
import argparse

# Locate the library file.
LIBFILE = 'core/libdisttbfast.so'
script_dir = os.path.dirname(os.path.realpath(__file__))
library_path = os.path.join(script_dir, LIBFILE)
if not os.path.isfile(library_path):
  ioe = IOError('Library file "'+LIBFILE+'" not found.')
  ioe.errno = errno.ENOENT
  raise ioe

disttbfast = ctypes.cdll.LoadLibrary(library_path)
disttbfast.disttbfast.restype = ctypes.c_int


def make_argparser():
  parser = argparse.ArgumentParser(description='Align a set of sequences.')
  parser.add_argument('input', type=argparse.FileType('r'), default=sys.stdin, nargs='?',
    help='Input sequences.')
  return parser


def main(argv):
  parser = make_argparser()
  args = parser.parse_args(argv[1:])
  seqs = []
  names = []
  name = None
  seq_num = 1
  for line_raw in args.input:
    line = line_raw.rstrip('\r\n')
    if line.startswith('>'):
      name = line[1:]
      continue
    else:
      seq = line
      if not name:
        name = str(seq_num)
      seqs.append(seq)
      names.append(name)
      seq_num += 1
      name = None
  alignment = align(seqs, names)
  print
  for i, seq in enumerate(alignment):
    print '>seq{}\n{}'.format(i, seq)


def align(seqs, names):
  ngui = len(seqs)
  assert ngui == len(names), (ngui, len(names))
  lgui = 0
  for seq in seqs:
    if len(seq) > lgui:
      lgui = len(seq)
  lgui *= 2
  seqgui = strlist_to_c(seqs)
  namegui = strlist_to_c(names)
  arglist = ('disttbfast', '-q', '0', '-E', '2', '-V', '-1.53', '-s', '0.0', '-W', '6', '-O', '-C',
             '0', '-D', '-b', '62', '-f', '-1.53', '-Q', '100.0', '-h', '0', '-F', '-X', '0.1')
  argv = strlist_to_c(arglist)
  argc = len(argv)
  disttbfast.disttbfast(ngui, lgui, namegui, seqgui, argc, argv, None)
  return seqgui


def strlist_to_c(strlist):
  c_strs = (ctypes.c_char_p * len(strlist))()
  for i, s in enumerate(strlist):
    c_strs[i] = ctypes.c_char_p(s)
  return c_strs


if __name__ == '__main__':
  sys.exit(main(sys.argv))
