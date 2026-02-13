%Skapa matriserna, gör så att storleken anapssas för antal noter
%Gör så att kan bestämma vilka sidor som har värme och hur hög
%For loop tills alla värden i matrisen är "ungefär samma"

N = 3; %roten ur totala antalet noder
n=1
M=zeros(N^2,N^2)
while n <= N^2

if mod(n, N) ~= 0
M(n,n+1)=0.25;

end

if mod(n-1, N) ~= 0
M(n,n-1)=0.25;

end

if n > N
M(n,n-N)=0.25;

end

if n < N^2-N %Fixa till för om finns noder på raden över)
M(n,n+N)=0.25;

end

n=n+1

end
M