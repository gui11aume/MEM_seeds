#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

'''This script is used to print the relatively complex transfer matrix
of the model with a main target and a duplicate.'''

class MatrixAsString:

   def __init__(self, gamma, k, p, kappa, N):
      self.gamma = gamma
      self.k     = k    
      self.p     = p
      self.kappa = kappa
      self.N     = N

   def power(self, txt, k):
      if k == 0: return '1'
      if k == 1: return txt
      return '(%s)^%d' % (txt, k)

   def AB_type_terms(self, k):
      terms = [self.power('%f*z' % (1-self.p), m) for m in range(k)]
      return '+'.join(terms)

   def tilde_AB_type_terms(self, k):
      qz = '%f*z' % (1-self.p)
      terms = [self.Cm(m) + '*' + self.power(qz, m) for m in range(k)]
      return '+'.join(terms)

   def A(self, k):
      cst_term = '%f*z*(1-(1-%f/3)^%d)' % (self.p, self.kappa, self.N)
      if k == 1: return cst_term
      return cst_term + '*(' + self.AB_type_terms(k) + ')'

   def B(self, k):
      cst_term = '%f*z*(1-%f/3)^%d' % (self.p, self.kappa, self.N)
      if k == 1: return cst_term
      return cst_term + '*(' + self.AB_type_terms(k) + ')'

   def Cm(self, k):
      if k == 0: return '1'
      if k == 1: return '(1-%f^%d)' % (self.kappa, self.N)
      return '(1-(1-(1-%f)^%d)^%d)' % (self.kappa, k, self.N)

   def tildeA(self, k):
      cst_term = '%f*z*(1-(1-%f/3)^%d)' % (self.p, self.kappa, self.N)
      if k == 1: return cst_term
      return cst_term + '*(' + self.tilde_AB_type_terms(k) + ')'

   def tildeB(self, k):
      cst_term = '%f*z*(1-%f/3)^%d' % (self.p, self.kappa, self.N)
      if k == 1: return cst_term
      return cst_term + '*(' + self.tilde_AB_type_terms(k) + ')'

   def r(self, k):
      if k == 1:
         return '(%s-%s)*%f*z' % (self.Cm(k), self.Cm(k-1), 1-self.p)
      return '(%s-%s)*(%f*z)^%d' % (self.Cm(k), self.Cm(k-1), 1-self.p, k)

   def show(self):
      rows = list()
      row1 = [self.tildeA(self.k), self.tildeB(self.k)] + \
                  [self.r(m) for m in range(1,self.gamma)]
      rows.append('[' + ','.join(row1) + ']')

      row2 = [self.tildeA(self.k), self.tildeB(self.gamma)] + \
                  [self.r(m) for m in range(1,self.gamma)]
      rows.append('[' + ','.join(row2) + ']')

      for i in range(self.gamma-1, 0, -1):
         rowi = [self.A(i), self.B(i)] + ['0'] * (self.gamma-1)
         rows.append('[' + ','.join(rowi) + ']')

      return '[' + ','.join(rows) + ']'


      

if __name__ == '__main__':
   M = MatrixAsString(3, 3, .01, .05, 1)
   print 'M := Matrix(' + M.show() + ')'
