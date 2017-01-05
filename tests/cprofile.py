
import cProfile
import unittest
import pstats

if __name__ == '__main__':
   suite = unittest.TestLoader().discover('.')
   def runtests():
      unittest.TextTestRunner().run(suite)
   s = cProfile.run('runtests()',sort='cumtime')
