%Skapa matriserna, gör så att storleken anapssas för antal noter
%Gör så att kan bestämma vilka sidor som har värme och hur hög
%For loop tills alla värden i matrisen är "ungefär samma"

n = 4 

classdef BasicClass
   properties
      Value {mustBeNumeric}
   end
   methods
      function r = roundOff(obj)
         r = round([obj.Value],2);
      end
      function r = multiplyBy(obj,n)
         r = [obj.Value]*n;
      end
   end
end

for n=4  %Antal noder



end

M=1/n.*[t1]