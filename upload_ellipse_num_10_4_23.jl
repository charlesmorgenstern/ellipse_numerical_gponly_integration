using Interpolations
using LaTeXStrings
using Roots
using Plots
using LinearAlgebra: norm, eigen, dot, cross
using Printf
using QuadGK: quadgk
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getpaths(n,a,h)
#generate gradient paths

b=sqrt(.5)*a

#allocate array of vectors for paths
gps=[Float64[] for a in 1:2*n]

#allocate array of start parameters
start=Matrix{Float64}(undef,n,1)
start[1]=0
start[end]=pi/2


r(t)=sqrt(a^2*sin(t)^2+b^2*cos(t)^2)
length, error = quadgk(r, 0, pi/2)
interval=length/(n-1)
tstep=pi/100000

for i=2:n-1
l=interval*(i-1) #get starting points
t=tstep
check=0
 while check<l
     check, error = quadgk(r,0,t)
     start[i]=t
     t+=tstep
 end
end
x=a*cos.(start)
y=b*sin.(start)


for i=1:n #RK4 for path generation
gpx=Float64[]
gpy=Float64[]
push!(gpx,x[i])
push!(gpy,y[i])
 while gpx[end]^2+2*gpy[end]^2<16
 yn=gpy[end]
 xn=gpx[end]
 ky1=4*yn
 ky2=4*(yn+h*ky1/2)
 ky3=4*(yn+h*ky2/2)
 ky4=4*(yn+h*ky3)
 yy=yn+h/6*(ky1+2*ky2+2*ky3+ky4)
 kx1=2*xn
 kx2=2*(xn+h*kx1/2)
 kx3=2*(xn+h*kx2/2)
 kx4=2*(xn+h*kx3)
 xx=xn+h/6*(kx1+2*kx2+2*kx3+kx4)
 push!(gpx,xx)
 push!(gpy,yy)
 end
y1=gpy[end-1]
y2=gpy[end]
x1=gpx[end-1]
x2=gpx[end]
m=(y2-y1)/(x2-x1) #force end points to be on ellipse
xx=(sqrt(2)*sqrt(-x1^2*m^2+2*x1*y1*m-y1^2+16*m^2+8)+2*x1*m^2-2*y1*m)/(2*m^2+1)
yy=m*(xx-x1)+y1
gpx[end]=xx
gpy[end]=yy
gps[2*i-1]=gpx
gps[2*i]=gpy
end
return gps
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function plotpaths(n,a,h)
#plot contours and gradient paths

f(x,y) = 1 - x^2 - 2*y^2 #the function
gps=getpaths(n,a,h)

# Create a grid of x and y values for plotting
x = range(0, 4.5, length=100)
y = range(0, 3.5, length=100)

# Evaluate the function on the grid
#z = [f(i, j) for i in x, j in y]
z = @. f(x', y)

# Create the contour plot
plt=contour(x, y, z, levels=40, color=:plasma, xlabel="x", ylabel="y")

t=collect(range(0,pi/2,200))
x=a*cos.(t)
y=a*sqrt(.5)*sin.(t)

plot!(x,y,legend=false)

x=4*cos.(t)
y=4*sqrt(.5)*sin.(t)

plot!(x,y,legend=false)

for i=1:n
x=gps[2*i-1]
y=gps[2*i]
plot!(x,y,legend=false)
end

return plt
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function curv(r)
#get curvature of level set at point r with exact formula
x=r[1]
y=r[2]


k=(2*x^2+4*y^2)/(x^2+4*y^2)^(3/2)

    return k
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function curvfdm2(r,h)
#get curvature of level set at point r with finite differences

f(x,y)=x^2+2*y^2 #function with ellipses for level sets

x2=r[1]
y2=r[2]
x1=x2-h
y1=y2-h
x3=x2+h
y3=y2+h
zx1=f(x1,y2)
zy1=f(x2,y1)
z2=f(x2,y2)
zx3=f(x3,y2)
zy3=f(x2,y3)

dx=(zx3-zx1)/(2*h)
dy=(zy3-zy1)/(2*h)
ddx=(zx3-2*z2+zx1)/(h^2)
ddy=(zy3-2*z2+zy1)/(h^2)
dxy=(f(x2+h,y2+h)-f(x2-h,y2+h)-f(x2+h,y2-h)+f(x2-h,y2-h))/(4*h^2)


k=(ddx*dy^2-2*dx*dxy*dy+ddy*dx^2)/(dx^2+dy^2)^(3/2)

end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getpathlengths(gps)
#calculate length of GP segments with Euclidean distance

n=floor(Int,(length(gps)/2))
lengths=[Float64[] for a in 1:n]

for j=1:n
x=gps[2*j-1]
y=gps[2*j]
nx=length(x)
h=Float64[]
for i=1:nx-1
push!(h,sqrt((x[i+1]-x[i])^2+(y[i+1]-y[i])^2))
end
lengths[j]=h
end
return lengths

end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getpathlengthsex(gps)
#caculate length of GP segments with integral

n=floor(Int,(length(gps)/2))
lengths=[Float64[] for a in 1:n]

for j=1:n-1
x=gps[2*j-1]
y=gps[2*j]
nx=length(x)
h=Float64[]
c=y[1]/x[1]^2
f(t)=sqrt(1+4*c^2*t^2)
for i=1:nx-1
len, error = quadgk(f, x[i], x[i+1])
push!(h,len)
end
lengths[j]=h
end

h=Float64[]
y=gps[2*n]
nx=length(y)
for i=1:nx-1
len=y[i+1]-y[i]
push!(h,len)
end
lengths[n]=h

return lengths

end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function checkpathlengths(n,a,h)
#calculate maximum relative error of GP segment lengths

gps=getpaths(n,a,h)
l1=getpathlengths(gps)
l2=getpathlengthsex(gps)
n=length(l1)
er=Float64[]

for i=1:n
temp=abs.(l1[i]-l2[i])/l2[i]
temp=maximum(maximum(temp))
push!(er,temp)
end

er=maximum(er)

return er

end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getcurvs(gps)
#calculate curvature of level sets at needed points with finite difference curv code
n=floor(Int,length(gps)/2)
curvs=[Float64[] for a in 1:n]
h=.001

for j=1:n
x=gps[2*j-1]
y=gps[2*j]
nx=length(x)
c=Float64[]
for i=1:nx
push!(c,curvfdm2((x[i],y[i]),h))
end
curvs[j]=c
end
return curvs
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getcurvsex(gps)
#calculate curvature of level sets at needed points with exact curvature formula
n=floor(Int,length(gps)/2)
curvs=[Float64[] for a in 1:n]


for j=1:n
x=gps[2*j-1]
y=gps[2*j]
scatter!(x,y)
nx=length(x)
c=Float64[]
for i=1:nx
x1=x[i]
y1=y[i]
if j==1
t=0
elseif j==n
t=pi/2
else
t=atan(y1/(sqrt(.5)*x1))
end
if j<n
a=x1/cos(t)
b=sqrt(.5)*a
else
a=y1/(sqrt(.5)*sin(t))
b=sqrt(.5)*a
end

#formula for curvature
temp=(a*b)/(sqrt(a^2*sin(t)^2+b^2*cos(t)^2))^3


push!(c,temp)
end

curvs[j]=c
end
return curvs

end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getcurvderivatives(gps)
#get derivative of level set curvature in direction of level set with finite difference
normal(x,y)=[-4*y,2*x]
n=floor(Int,length(gps)/2)
dk=[Float64[] for a in 1:n]
h=.00001

for i=1:n
d=Float64[]
n2=length(gps[2*i])
for j=1:n2
x1=gps[2*i-1][j]
y1=gps[2*i][j]
kx1=curv((x1-h,y1))
kx2=curv((x1+h,y1))
ky1=curv((x1,y1-h))
ky2=curv((x1,y1+h))
dx=(kx2-kx1)/(2*h)
dy=(ky2-ky1)/(2*h)
n1=normal(x1,y1)
n1=n1./sqrt(n1[1]^2+n1[2]^2)
dk2=n1[1]*dx+n1[2]*dy

push!(d,dk2)
end
dk[i]=d
end

return dk
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getcurvderivativesex(gps)
#get derivative of level set curvature in direction of level set with exact formula
nu(x,y)=[-4*y,2*x]
n=floor(Int,length(gps)/2)
dk=[Float64[] for a in 1:n]
h=.0000001
dx(x,y)=-(2*x*(x^2-2*y^2))/(x^2+4*y^2)^(5/2)
dy(x,y)=-(16*y*(x^2+y^2))/(x^2+4*y^2)^(5/2)

for i=1:n
d=Float64[]
n2=length(gps[2*i])
for j=1:n2
x1=gps[2*i-1][j]
y1=gps[2*i][j]
dx1=dx(x1,y1)
dy1=dy(x1,y1)
n1=nu(x1,y1)
n1=n1./sqrt(n1[1]^2+n1[2]^2)
dk2=n1[1]*dx1+n1[2]*dy1

push!(d,dk2)
end
dk[i]=d
end

return dk
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getpathcurv(gps)
#get curvature of GP with finite differences
n=floor(Int,length(gps)/2)
k=[Float64[] for a in 1:n-1]
#erk=Float64[]
h2=(gps[1][2]-gps[1][1])/(2*gps[1][1])


for i=1:n-1
np=length(gps[2*i-1])
k2=Float64[]
for j=1:np
 x2=gps[2*i-1][j]
 y2=gps[2*i][j]
   ky1=4*y2
   ky2=4*(y2+h2*ky1/2)
   ky3=4*(y2+h2*ky2/2)
   ky4=4*(y2+h2*ky3)
   y3=y2+h2/6*(ky1+2*ky2+2*ky3+ky4)
   kx1=2*x2
   kx2=2*(x2+h2*kx1/2)
   kx3=2*(x2+h2*kx2/2)
   kx4=2*(x2+h2*kx3)
   x3=x2+h2/6*(kx1+2*kx2+2*kx3+kx4)
   ky1=-4*y2
   ky2=-4*(y2+h2*ky1/2)
   ky3=-4*(y2+h2*ky2/2)
   ky4=-4*(y2+h2*ky3)
   y1=y2+h2/6*(ky1+2*ky2+2*ky3+ky4)
   kx1=-2*x2
   kx2=-2*(x2+h2*kx1/2)
   kx3=-2*(x2+h2*kx2/2)
   kx4=-2*(x2+h2*kx3)
   x1=x2+h2/6*(kx1+2*kx2+2*kx3+kx4)

 dx=(x3-x1)/2
 dy=(y3-y1)/2
 ddx=x3-2*x2+x1
 ddy=y3-2*y2+y1
temp=(dx*ddy-dy*ddx)/(dx^2+dy^2)^(3/2)
push!(k2,temp)
end
k[i]=k2
end

return k
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getpathcurvex(gps)
#get curvature of GP with exact formula
n=floor(Int,length(gps)/2)
k=[Float64[] for a in 1:n-1]
#erk=Float64[]
h2=(gps[1][2]-gps[1][1])/(2*gps[1][1])


for i=1:n-1
np=length(gps[2*i-1])
k2=Float64[]
for j=1:np
 x=gps[2*i-1][j]
 y=gps[2*i][j]
if i==1
c=0
else
 c=y/x^2
end
temp=(2*c)/((1+4*c^2*x^2)^(3/2)) 
push!(k2,temp)
end
k[i]=k2
end
return k
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getleftpathcurv(gps)
#get GP curvature of neighboring GP at corresponding points
n=floor(Int,length(gps)/2)
k=[Float64[] for a in 1:n-1]
gps2=[Float64[] for a in 1:2*n-2]
h=(gps[1][2]-gps[1][1])/(2*gps[1][1])
h2=h
h=h/10000

f(x,y) = 1.0 -x^2-2*y^2

for i=1:n-1

xs=Float64[]
ys=Float64[]
np=length(gps[2*i-1])
push!(xs,gps[2*i+1][1])
push!(ys,gps[2*i+2][1])


for j=2:np-1
 x=gps[2*i-1][j]
 y=gps[2*i][j]
 val=f(x,y)
 m=2
 x2=gps[2*i+1][2]
 y2=gps[2*i+2][2]
 check=f(x2,y2)
  while check>val
   m=m+1
   x2=gps[2*i+1][m]
   y2=gps[2*i+2][m]
   check=f(x2,y2)
  end
 m=m-1
 x2=gps[2*i+1][m]
 y2=gps[2*i+2][m]
 check=f(x2,y2)
  while check>val
   ky1=4*y2
   ky2=4*(y2+h*ky1/2)
   ky3=4*(y2+h*ky2/2)
   ky4=4*(y2+h*ky3)
   y2=y2+h/6*(ky1+2*ky2+2*ky3+ky4)
   kx1=2*x2
   kx2=2*(x2+h*kx1/2)
   kx3=2*(x2+h*kx2/2)
   kx4=2*(x2+h*kx3)
   x2=x2+h/6*(kx1+2*kx2+2*kx3+kx4)
   check=f(x2,y2)
  end
 push!(xs,x2)
 push!(ys,y2)
end
gps2[2*i-1]=xs
gps2[2*i]=ys
end

for i=1:n-1
k2=Float64[]
np=length(gps2[2*i-1])
for j=1:np
 x2=gps2[2*i-1][j]
 y2=gps2[2*i][j]
   ky1=4*y2
   ky2=4*(y2+h2*ky1/2)
   ky3=4*(y2+h2*ky2/2)
   ky4=4*(y2+h2*ky3)
   y3=y2+h2/6*(ky1+2*ky2+2*ky3+ky4)
   kx1=2*x2
   kx2=2*(x2+h2*kx1/2)
   kx3=2*(x2+h2*kx2/2)
   kx4=2*(x2+h2*kx3)
   x3=x2+h2/6*(kx1+2*kx2+2*kx3+kx4)
   ky1=-4*y2
   ky2=-4*(y2+h2*ky1/2)
   ky3=-4*(y2+h2*ky2/2)
   ky4=-4*(y2+h2*ky3)
   y1=y2+h2/6*(ky1+2*ky2+2*ky3+ky4)
   kx1=-2*x2
   kx2=-2*(x2+h2*kx1/2)
   kx3=-2*(x2+h2*kx2/2)
   kx4=-2*(x2+h2*kx3)
   x1=x2+h2/6*(kx1+2*kx2+2*kx3+kx4)

 dx=(x3-x1)/2
 dy=(y3-y1)/2
 ddx=x3-2*x2+x1
 ddy=y3-2*y2+y1
temp=(dx*ddy-dy*ddx)/(dx^2+dy^2)^(3/2)
push!(k2,temp)

end
k[i]=k2

end

return k

end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getleftpathlengths(gps)
#get lengths of corresponding GP segments of neighboring GP
n=floor(Int,length(gps)/2)
k=[Float64[] for a in 1:n-1]
gps2=[Float64[] for a in 1:2*n-2]
h=(gps[1][2]-gps[1][1])/(2*gps[1][1])
h2=h
h=h/10000

f(x,y) = 1.0 -x^2-2*y^2

for i=1:n-1

xs=Float64[]
ys=Float64[]
np=length(gps[2*i-1])
push!(xs,gps[2*i+1][1])
push!(ys,gps[2*i+2][1])


for j=2:np-1
 x=gps[2*i-1][j]
 y=gps[2*i][j]
 val=f(x,y)
 m=2
 x2=gps[2*i+1][2]
 y2=gps[2*i+2][2]
 check=f(x2,y2)
  while check>val
   m=m+1
   x2=gps[2*i+1][m]
   y2=gps[2*i+2][m]
   check=f(x2,y2)
  end
 m=m-1
 x2=gps[2*i+1][m]
 y2=gps[2*i+2][m]
 check=f(x2,y2)
  while check>val
   ky1=4*y2
   ky2=4*(y2+h*ky1/2)
   ky3=4*(y2+h*ky2/2)
   ky4=4*(y2+h*ky3)
   y2=y2+h/6*(ky1+2*ky2+2*ky3+ky4)
   kx1=2*x2
   kx2=2*(x2+h*kx1/2)
   kx3=2*(x2+h*kx2/2)
   kx4=2*(x2+h*kx3)
   x2=x2+h/6*(kx1+2*kx2+2*kx3+kx4)
   check=f(x2,y2)
  end
 push!(xs,x2)
 push!(ys,y2)
end
push!(xs,gps[2*i+1][end])
push!(ys,gps[2*i+2][end])
gps2[2*i-1]=xs
gps2[2*i]=ys
end


n=floor(Int,(length(gps2)/2))
lengths=[Float64[] for a in 1:n]

for j=1:n
x=gps2[2*j-1]
y=gps2[2*j]
nx=length(x)
h=Float64[]
for i=1:nx-1
push!(h,sqrt((x[i+1]-x[i])^2+(y[i+1]-y[i])^2))
end
lengths[j]=h
end
return lengths

end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getarclengths(gps,curvs,lengths,lengths2,dk,a,k,k2)
#get arc length of level set segments
n=floor(Int,length(gps)/2)


#allocate array of vectors for arc lengths
arcs=[Float64[] for a in 1:n-1]

f(t)=sqrt(a^2*sin(t)^2+a^2*.5*cos(t)^2)

for i=1:n-1
aa=Float64[]
n2=length(gps[2*i])
x1=gps[2*i-1][1]
x2=gps[2*i+1][1]
t1=acos(x1/a)
t2=acos(x2/a)
l0, error = quadgk(f, t1, t2)
push!(aa,l0)
for j=2:n2 #recurrence relation
li=aa[j-1]*(1+curvs[i][j-1]*lengths[i][j-1]+dk[i][j-1]/2*lengths[i][j-1]*aa[j-1])-k[i][j-1]*lengths[i][j-1]^2/2+k2[i][j-1]*lengths2[i][j-1]^2/2
push!(aa,li)
end
arcs[i]=aa
end

return arcs

end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getarea(n,a,h)
#calculate area of ellipse

gps=getpaths(n,a,h)
curvs=getcurvs(gps)
lengths=getpathlengths(gps)
lengths2=getleftpathlengths(gps)
dk=getcurvderivatives(gps)
k=getpathcurv(gps)
k2=getleftpathcurv(gps)
arcs=getarclengths(gps,curvs,lengths,lengths2,dk,a,k,k2)

area=0.0
for i=1:n-1
n2=length(gps[2*i])
for j=1:n2-1
area+=lengths[i][j]*(arcs[i][j]) #left hand
end
end

return area
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getareagb(n,a,h)
#calculate area of individual GBs

gps=getpaths(n,a,h)
curvs=getcurvs(gps)
lengths=getpathlengths(gps)
lengths2=getleftpathlengths(gps)
dk=getcurvderivatives(gps)
k=getpathcurv(gps)
k2=getleftpathcurv(gps)
arcs=getarclengths(gps,curvs,lengths,lengths2,dk,a,k,k2)
areas=Matrix{Float64}(undef,n-1,1)

for i=1:n-1
area=0.0
n2=length(gps[2*i])
for j=1:n2-1
area+=lengths[i][j]*(arcs[i][j]) #left hand
end
areas[i]=area
end

return areas
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getareagreens(n,a,h)
#calculate area of individual GBs with Green's Thm

gps=getpaths(n,a,h)
areas=Matrix{Float64}(undef,n-1,1)

for i=1:n-1

px1=gps[2*i-1]
py1=gps[2*i]
px3=gps[2*i+1]
py3=gps[2*i+2]

aa=4
b=aa*sqrt(.5)

if i==1
t1=0
else
t1=acos(px1[end]/aa)
end
t2=acos(px3[end]/aa)
t=collect(range(t1,t2,1000))

px2=aa*cos.(t)
py2=b*sin.(t)

b=a*sqrt(.5)
if i==1
t1=0
else
t1=acos(px1[1]/a)
end
t2=acos(px3[1]/a)
t=collect(range(t1,t2,1000))

px4=a*cos.(t)
py4=b*sin.(t)

px4=reverse(px4,dims=1)
py4=reverse(py4,dims=1)
px3=reverse(px3,dims=1)
py3=reverse(py3,dims=1)

n1=length(px1)
n2=length(px2)
n3=length(px3)
n4=length(px4)

i1=0.0
for j=2:n1
x=px1[j]
y=py1[j]
x1=px1[j-1]
y1=py1[j-1]
h=sqrt((x-x1)^2+(y-y1)^2)
dy=(y-y1)/h
dx=(x-x1)/h
i1+=(x*dy-y*dx)*h
end

i2=0.0
for j=2:n2
x=px2[j]
y=py2[j]
x1=px2[j-1]
y1=py2[j-1]
h=sqrt((x-x1)^2+(y-y1)^2)
dy=(y-y1)/h
dx=(x-x1)/h
i2+=(x*dy-y*dx)*h
end

i3=0.0
for j=2:n3
x=px3[j]
y=py3[j]
x1=px3[j-1]
y1=py3[j-1]
h=sqrt((x-x1)^2+(y-y1)^2)
dy=(y-y1)/h
dx=(x-x1)/h
i3+=(x*dy-y*dx)*h
end

i4=0.0
for j=2:n4
x=px4[j]
y=py4[j]
x1=px4[j-1]
y1=py4[j-1]
h=sqrt((x-x1)^2+(y-y1)^2)
dy=(y-y1)/h
dx=(x-x1)/h
i4+=(x*dy-y*dx)*h
end

area=.5*(i1+i2+i3+i4)
areas[i]=area
end

return areas
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getareaerror(n,a,h)
#get relative error of area of ellipse

area=getarea(n,a,h)
b=sqrt(.5)*a
exact=4*sqrt(.5)*pi-a*b*pi/4

err=abs(exact-area)/exact
display("exact")
display(exact)
display("area")
display(area)

return err
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getn(gps)
#get maximum number of points along a GP 

ngp=floor(Int,length(gps)/2)

n=Float64[]

for i=1:ngp
len=length(gps[2*i])
push!(n,len)
end

n=maximum(n)

return n
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function errortablearea(a,nn,nstart,hstart)
#create relative error table for area of ellipse with increasing number of GPs

n=Matrix{Int64}(undef,nn,1)
h=Matrix{Float64}(undef,nn,1)
n[1]=nstart
h[1]=hstart

for i=2:nn
n[i]=n[i-1]*2
h[i]=h[i-1]/2
end

err=Matrix{Float64}(undef,nn,1)
eoc=Matrix{Float64}(undef,nn,1)
narcs=Matrix{Float64}(undef,nn,1)

gps=getpaths(n[1],a,h[1])
narcs[1]=getn(gps)
err[1]=getareaerror(n[1],a,h[1])
for i=2:nn
display("working on:")
display(i)
display("ngps:")
display(n[i])
display("h:")
display(h[i])
gps=getpaths(n[i],a,h[i])
narcs[i]=getn(gps)
err[i]=getareaerror(n[i],a,h[i])
eoc[i]=log(err[i-1]/err[i])/log(2)
end

@printf "\n Relative error of area using left hand"
@printf "\n ngps    narcs      rel. error                eoc"
@printf "\n %g       %g         %g                 n/a" n[1] narcs[1]  err[1]
for i=2:nn
@printf "\n %g       %g            %g              %g" n[i] narcs[i] err[i] eoc[i]
end


end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function plotarclengtherror(n,a,h)
#plot relative error of arc lengths on largest ellipse

gps=getpaths(n,a,h)
curvs=getcurvs(gps)
lengths=getpathlengths(gps)
lengths2=getleftpathlengths(gps)
dk=getcurvderivatives(gps)
k=getpathcurv(gps)
k2=getleftpathcurv(gps)
arcs=getarclengths(gps,curvs,lengths,lengths2,dk,a,k,k2)

f(t)=sqrt(16*sin(t)^2+8*cos(t)^2)
err=Matrix{Float64}(undef,n-1,1)
t=Matrix{Float64}(undef,n-1,1)

for i=1:n-1
x1=gps[2*i-1][end]
x2=gps[2*i+1][end]
if x1>4
x1=4
end
t1=acos(x1/4)
t2=acos(x2/4)
t[i]=(t1+t2)/2
exact, error = quadgk(f,t1,t2)
length=arcs[i][end]
err[i]=abs(exact-length)/exact
end
fig=plot(t,err,legend=false,xlabel="Center Angle of Arc in Radians",ylabel="Relative Error of Arc Length",title="x-intercept of first ellipse=$a",ylims=(0,.1))

return fig
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function animatearclengtherror()
#animate plots of arc length error as starting ellipse size changes
na=5
a=collect(range(.5,3.9,na))
n=640
h=.00078125
anim = @animate for i ∈ 1:na
plotarclengtherror(n,a[i],h)
end
gif(anim, "arcerror_fps5.gif", fps = 5)
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function animatesetup()
#animate gradient paths as starting ellipse size changes

na=7
a=[.5,1,1.5,2,2.5,3,3.5]
n=80
h=.001
anim = @animate for i ∈ 1:na
plotpaths(n,a[i],h)
end
gif(anim, "setup_fps5.gif", fps = 3)
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getc(gps)
#get parameters c of parabolas y=cx^2 
ngp=floor(Int,length(gps)/2)
c=Matrix{Float64}(undef,ngp,1)

for i=1:ngp
x=gps[2*i-1][1]
y=gps[2*i][1]
c[i]=y/x^2
end

return c
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getgperrorL2(gps)
#get L2 error of gradient paths to check RK4

lengths=getpathlengths(gps)
c=getc(gps)
ngp=floor(Int,length(gps)/2)
L2=Matrix{Float64}(undef,ngp,1)

for i=1:ngp
n=length(gps[2*i-1])
L2[i]=0

for j=1:n-1
x=gps[2*i-1][j]
y=gps[2*i][j]
a=sqrt(x^2+2*y^2)
if i==1
yex=0.0
else
 f(y)=2*y^2+y/c[i]-a^2
 yex=fzero(f,y)
end
if i==1
xex=a
else
xex=sqrt(yex/c[i])
end
L2[i]=L2[i]+((x-xex)^2+(y-yex)^2)*lengths[i][j]
end
L2[i]=sqrt(L2[i])

end
return L2
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function plotpathswexact(n,a,h)
#plot RK4 GPs and exact GPs on same plot

f(x,y) = 1 - x^2 - 2*y^2 #the function
gps=getpaths(n,a,h)

# Create a grid of x and y values for plotting
x = range(0, 4.5, length=100)
y = range(0, 3.5, length=100)

# Evaluate the function on the grid
#z = [f(i, j) for i in x, j in y]
z = @. f(x', y)

# Create the contour plot
plt=contour(x, y, z, levels=40, color=:plasma, xlabel="x", ylabel="y")

t=collect(range(0,pi/2,200))
x=a*cos.(t)
y=a*sqrt(.5)*sin.(t)

plot!(x,y,legend=false)

x=4*cos.(t)
y=4*sqrt(.5)*sin.(t)

plot!(x,y,legend=false)

for i=1:n
x=gps[2*i-1]
y=gps[2*i]
plot!(x,y,legend=false)
end

c=getc(gps)
for i=1:n
x=collect(range(0,gps[2*i-1][end],100))
y=c[i].*x.^2
plot!(x,y,legend=false,c=:black)
end
return plt
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getarclengthsex(gps)
#get exact arc lengths of level sets with integral

c=getc(gps)
ngp=floor(Int,length(gps)/2)
arcs=[Float64[] for a in 1:ngp-1]

for i=1:ngp-1
n=length(gps[2*i-1])
aa=Float64[]
for j=1:n
x=gps[2*i-1][j]
y=gps[2*i][j]
a=sqrt(x^2+2*y^2)
t=acos(x/a)
 f(y)=2*y^2+y/c[i+1]-a^2
 y2=fzero(f,y)
 x2=sqrt(y2/c[i+1])
t2=acos(x2/a)
f2(t)=sqrt(a^2*sin(t)^2+.5*a^2*cos(t)^2)
exact, error = quadgk(f2,t,t2)
push!(aa,exact)
end
arcs[i]=aa
end

return arcs
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getarclengtherrors(n,a,h)
#get relative error of level set arc lengths

gps=getpaths(n,a,h)
curvs=getcurvs(gps)
lengths=getpathlengths(gps)
lengths2=getleftpathlengths(gps)
dk=getcurvderivatives(gps)
k=getpathcurv(gps)
k2=getleftpathcurv(gps)
arcs=getarclengths(gps,curvs,lengths,lengths2,dk,a,k,k2)
arcsex=getarclengthsex(gps)
ngp=floor(Int,length(gps)/2)
err=[Float64[] for a in 1:ngp-1]

for i=1:ngp-1
n=length(gps[2*i-1])
e=Float64[]
for j=1:n
er=abs(arcs[i][j]-arcsex[i][j])/arcsex[i][j]
push!(e,er)
end
err[i]=e
end
return err
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function plotarclengtherrors(n,a,h)
#plot the relative error of the arc lengths

err=getarclengtherrors(n,a,h)
f(x,y) = 1 - x^2 - 2*y^2 #the function
gps=getpaths(n,a,h)
# Create a grid of x and y values for plotting
x = range(0, 4.5, length=100)
y = range(0, 3.5, length=100)

# Evaluate the function on the grid
#z = [f(i, j) for i in x, j in y]
z = @. f(x', y)

# Create the contour plot
plt=contour(x, y, z, levels=40, color=:plasma, xlabel="x", ylabel="y")

t=collect(range(0,pi/2,200))
x=a*cos.(t)
y=a*sqrt(.5)*sin.(t)

plot!(x,y,legend=false)

x=4*cos.(t)
y=4*sqrt(.5)*sin.(t)

plot!(x,y,legend=false)

for i=1:n-1
x=gps[2*i-1]
y=gps[2*i]
e=250*err[i]
scatter!(x,y,markersize=e,legend=false)
end
display(maximum(maximum(err)))
return plt

end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function plotgberror(n,a,h)
#plot relative error of individual GBs area title with number of GBs
ex=getareagreens(n,a,h)
ap=getareagb(n,a,h)
er=abs.(ex-ap)./abs.(ex)
h=(pi/2)/(n-1)
t=collect(range(0+.5*h,pi/2-.5*h,n-1))
fig=plot(t,er,legend=false,xlabel="center angle of GB on starting ellipse",ylabel="relative error of area",ylims = (0,.03),title = L"Number  of  GBs = %$(n)")
return fig
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function animategberror()
#animate plots of relative error of individual GBs area with varying number of GBs
na=7
n=[10,20,40,80,160,320,640]
h=[.01,.01*.5^1,.01*.5^2,.01*.5^3,.01*.5^4,.01*.5^5,.01*.5^6]
a=3.5

anim = @animate for i ∈ 1:na
plotgberror(n[i],a,h[i])
end

gif(anim, "gberror3p5_fps1.gif", fps = 1)
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function plotgberror2(n,a,h)
#plot relative error of individual GBs area title with starting ellipse size
ex=getareagreens(n,a,h)
ap=getareagb(n,a,h)
er=abs.(ex-ap)./abs.(ex)
h=(pi/2)/(n-1)
t=collect(range(0+.5*h,pi/2-.5*h,n-1))
fig=plot(t,er,legend=false,xlabel="center angle of GB on starting ellipse",ylabel="relative error of area",ylims = (0,.015),title = L"x-intercept\; of\; starting\; ellipse = %$(a)")
return fig
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function animategberror2()
#animate plots of relative error of individual GBs area with varying starting ellipse size
na=7
n=160 
h=.01*.5^4 
a=[.5,1,1.5,2,2.5,3,3.5]

anim = @animate for i ∈ 1:na
plotgberror2(n,a[i],h)
end

gif(anim, "gberrora_fps1.gif", fps = 1)
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function plotcurverrors(n,a,h)
#plot errors of level set curvatures
err=getcurverrors(n,a,h)
f(x,y) = 1 - x^2 - 2*y^2 #the function
gps=getpaths(n,a,h)

# Create a grid of x and y values for plotting
x = range(0, 4.5, length=100)
y = range(0, 3.5, length=100)

# Evaluate the function on the grid
#z = [f(i, j) for i in x, j in y]
z = @. f(x', y)

# Create the contour plot
plt=contour(x, y, z, levels=40, color=:plasma, xlabel="x", ylabel="y")

t=collect(range(0,pi/2,200))
x=a*cos.(t)
y=a*sqrt(.5)*sin.(t)

plot!(x,y,legend=false)

x=4*cos.(t)
y=4*sqrt(.5)*sin.(t)

plot!(x,y,legend=false)

for i=1:n-1
x=gps[2*i-1]
y=gps[2*i]
e=2.5*10^9*err[i]
scatter!(x,y,markersize=e,legend=false)
end
display(maximum(maximum(err)))

return plt

end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getcurverrors(n,a,h)
#get relative error of level set curvatures
gps=getpaths(n,a,h)
curvs=getcurvs(gps)
curvsex=getcurvsex(gps)
ngp=floor(Int,length(gps)/2)
err=[Float64[] for a in 1:ngp-1]

for i=1:ngp-1
n=length(gps[2*i-1])
e=Float64[]
for j=1:n
er=abs(curvs[i][j]-curvsex[i][j])/abs(curvsex[i][j])
push!(e,er)
end
err[i]=e
end
return err
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getcurvsex(gps)
#get exact level set curvatures
n=floor(Int,length(gps)/2)
curvs=[Float64[] for a in 1:n]
h=.001

for j=1:n
x=gps[2*j-1]
y=gps[2*j]
nx=length(x)
c=Float64[]
for i=1:nx
push!(c,curv((x[i],y[i])))
end
curvs[j]=c
end
return curvs
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getpatherrors(n,a,h)
#get relative error of GP segment lengths
gps=getpaths(n,a,h)
len=getpathlengths(gps)
lenex=getpathlengthsex(gps)
ngp=floor(Int,length(gps)/2)
err=[Float64[] for a in 1:ngp-1]

for i=1:ngp-1
n=length(gps[2*i-1])
e=Float64[]
for j=1:n-1
er=abs(len[i][j]-lenex[i][j])/abs(lenex[i][j])
push!(e,er)
end
err[i]=e
end
return err
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function plotpatherrors(n,a,h)
#plot error of GP segment lengths
err=getpatherrors(n,a,h)
f(x,y) = 1 - x^2 - 2*y^2 #the function
gps=getpaths(n,a,h)

# Create a grid of x and y values for plotting
x = range(0, 4.5, length=100)
y = range(0, 3.5, length=100)

# Evaluate the function on the grid
#z = [f(i, j) for i in x, j in y]
z = @. f(x', y)

# Create the contour plot
plt=contour(x, y, z, levels=40, color=:plasma, xlabel="x", ylabel="y")

t=collect(range(0,pi/2,200))
x=a*cos.(t)
y=a*sqrt(.5)*sin.(t)

plot!(x,y,legend=false)

x=4*cos.(t)
y=4*sqrt(.5)*sin.(t)

plot!(x,y,legend=false)
for i=1:n-1
x=gps[2*i-1]
y=gps[2*i]
e=5*10^5*err[i][1:end-1]
scatter!(x,y,markersize=e,legend=false)
end
display(maximum(maximum(err)))

return plt

end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getpathcurverrors(n,a,h)
#get relative error of GP curvatures
gps=getpaths(n,a,h)
k=getpathcurv(gps)
kex=getpathcurvex(gps)
ngp=floor(Int,length(gps)/2)
err=[Float64[] for a in 1:ngp-1]

for i=2:ngp-1
n=length(gps[2*i-1])
e=Float64[]
for j=1:n-1
er=abs(k[i][j]-kex[i][j])/abs(kex[i][j])
push!(e,er)
end
err[i]=e
end
return err
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function plotpathcurverrors(n,a,h)
#plot error of GP curvatures
err=getpathcurverrors(n,a,h)
f(x,y) = 1 - x^2 - 2*y^2 #the function
gps=getpaths(n,a,h)

# Create a grid of x and y values for plotting
x = range(0, 4.5, length=100)
y = range(0, 3.5, length=100)

# Evaluate the function on the grid
#z = [f(i, j) for i in x, j in y]
z = @. f(x', y)

# Create the contour plot
plt=contour(x, y, z, levels=40, color=:plasma, xlabel="x", ylabel="y")

t=collect(range(0,pi/2,200))
x=a*cos.(t)
y=a*sqrt(.5)*sin.(t)

plot!(x,y,legend=false)

x=4*cos.(t)
y=4*sqrt(.5)*sin.(t)

plot!(x,y,legend=false)
for i=2:n-1
x=gps[2*i-1]
y=gps[2*i]
e=.25*10^4*err[i]
scatter!(x,y,markersize=e,legend=false)
end
display(maximum(maximum(err[2:end])))

return plt

end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function getcurvderivativeerrors(n,a,h)
#get relative error of level set curvature derivatives
gps=getpaths(n,a,h)
dk=getcurvderivatives(gps)
dkex=getcurvderivativesex(gps)
ngp=floor(Int,length(gps)/2)
err=[Float64[] for a in 1:ngp-1]

for i=2:ngp-1
n=length(gps[2*i-1])
e=Float64[]
for j=1:n
er=abs(dk[i][j]-dkex[i][j])/abs(dkex[i][j])
push!(e,er)
end
err[i]=e
end
return err
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
function plotcurvderivativeerrors(n,a,h)
#plot error of level set curvature derivatives
err=getcurvderivativeerrors(n,a,h)
f(x,y) = 1 - x^2 - 2*y^2 #the function
gps=getpaths(n,a,h)

# Create a grid of x and y values for plotting
x = range(0, 4.5, length=100)
y = range(0, 3.5, length=100)

# Evaluate the function on the grid
#z = [f(i, j) for i in x, j in y]
z = @. f(x', y)

# Create the contour plot
plt=contour(x, y, z, levels=40, color=:plasma, xlabel="x", ylabel="y")

t=collect(range(0,pi/2,200))
x=a*cos.(t)
y=a*sqrt(.5)*sin.(t)

plot!(x,y,legend=false)

x=4*cos.(t)
y=4*sqrt(.5)*sin.(t)

plot!(x,y,legend=false)
for i=2:n-2
x=gps[2*i-1]
y=gps[2*i]
e=10^9*err[i][1:end]
scatter!(x,y,markersize=e,legend=false)
end
display(maximum(maximum(err)))

return plt

end
