import argparse

parser = argparse.ArgumentParser(prog='PROG')


subparsers = parser.add_subparsers(help='sub-command help')

parser_main = subparsers.add_parser('', help='a help')
parser_a.add_argument('bar', type=int, help='bar help')

# create the parser for the "b" command
parser_b = subparsers.add_parser('b', help='b help')
parser_b.add_argument('--baz', choices='XYZ', help='baz help')