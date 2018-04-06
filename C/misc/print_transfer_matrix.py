#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

'''This script is used to print the relatively complex transfer matrix
of the model with a main target and a duplicate.'''

class MatrixAsString:

   def sqbr(txt):
      return '[' + txt + ']'

   def rdbr(txt):
      return '(' + txt + ')'

   def pow(txt, m):
      if m == 0: return '1'
      if m == 1: return txt
      return self.rdbr('%s^%d' % txt, m)

   def A(k, p, kappa, N):
      cst_term = '%f*(1-(1-%f/3)^%d)*z' % (p, kappa, N)
      if k == 1: return cst_term

      other_terms = [self.pow('%f*z' % 1-p, m) for m in range(k)]
      self.sqbr('+'.join(other_terms))


   def first_term(n):
      '''Return a string of the first term of the row.'''
      if n == 0: return 'd*z'
      return 'd*z*(' + '+'.join([az(a) for a in range(n+1)]) + ')'

   def shift_right(L, shift):
      '''Shift the terms to the right.'''
      if shift == 0: return L
      return ['0']*shift + L[:-shift]

   def overwrite_tail(L, n):
      '''Replace the tail elements by 0.'''
      if n == 0: return L
      return L[:-n] + ['0']*n 

   # Two complete lists of terms.
   L1 = ['c*z'] + ['c*z*%s' % az(a) for a in range(1,d-1)]
   L2 = ['b*z'] + ['b*z*%s' % az(a) for a in range(1,d-1)]

   # First row.
   row_terms = [first_term(d-1)] + L1 + L2
   rows = [brackets(', '.join(row_terms))]

   # Next series of rows.
   for i in range(1,d):
      row_terms = [first_term(d-1-i)] + \
            shift_right(L1,i) + overwrite_tail(L2,i-1)
      rows.append(brackets(', '.join(row_terms)))

   # Next series of rows.
   for i in range(1,d):
      row_terms = [first_term(d-1-i)] + \
            overwrite_tail(L1,i-1) + shift_right(L2,i)
      rows.append(brackets(', '.join(row_terms)))

   return brackets(', '.join(rows))

def write_assignment_line(d):
   head = 'M%d := Matrix(' %d
   tail = '):'
   return head + create_matrix_as_string(d) + tail

if __name__ == '__main__':
   print write_assignment_line(int(sys.argv[1]))
