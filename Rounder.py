import math as math

class Rounder:
   'Common base class for all Samples'

   def getExponent(self, s):
  
       exponent = 0 
       v = math.fabs(s) 
       if(v > 1):
           while v > math.pow(10, exponent):
                exponent = exponent + 1
           exponent = exponent - 1 
       else:
           while v < math.pow(10, exponent):
                exponent = exponent - 1
                #v = v * 10.0
       return exponent  

 
   def getRound(self, d, expo):

       if(math.fabs(expo)>10):
         print "Is it really such a weird number?"
         return 0 
       val = math.floor(d/math.pow(10, expo))*math.pow(10, expo)
       val1 = (1+math.floor(d/math.pow(10, expo)))*math.pow(10, expo)
                
       if (2*d < val1+val):
           return val
       else:
           return val1

   def roundTo(self, a, sa):

       exposigma = self.getExponent(sa)
       sigma = self.getRound(sa, exposigma)
       val = self.getRound(a, exposigma)
       return [val, sigma]

   def toString(self, a):

       exposigma = self.getExponent(a)
       return str(self.getRound(a, exposigma))

   def toStringB(self, a, sa):

       p = self.roundTo(a, sa)
       return str(p[0]) + " +/- " + str(p[1])
