include("graph.jl")
include("edgecore.jl")
using LinearAlgebra
#using Laplacians


fname = open("filename.txt", "r")
fout = open("ans.txt", "w")
str   = readline(fname);
n     = parse(Int, str);
k=20;

    str = readline(fname);
    str = split(str);
    G   = get_graph(str[1]);
n=G.n;

L=lap(G);
selc1=zeros(k);
selc2=zeros(2k);
d=zeros(n);
for i=1:n
    d[i]=L[i,i];
end
x=argmax(d);
for i=1:k
    cho=argmax(d);
    selc1[i]=cho;
    d[cho]=-1;
end
f=eigen(L);

A=-L;
for i=1:n
    A[i,i]=0;
end
#f=eigen(A);
phi=zeros(n);
#sumlmd=zeros(n);
u=f.vectors;
for i=1:n
#    sumlmd[i]=sum(u[:,i]);
end
for i=1:n
    if u[x,i]<0
        for j=1:n
            #u[j,i]=-u[j,i];
        end
    end
end

for i=1:n
    for j=1:n
#        phi[i]+=(f.values[j])*(f.vectors[i,j])*f.vectors[j,j]/abs(f.vectors[j,j]);
        phi[i]+=(f.values[j]^2)*(u[i,j]);
    end
end
for i=1:n
#    phi[i]=abs(phi[i]);
end

for i=1:k
    cho=argmax(phi);
    selc2[i]=cho;
    phi[cho]=-10000000;
#    cho=argmin(phi);
#    selc2[i*2-1]=cho;
#    phi[cho]=0;
end


L=lap(G);
I=zeros(n,n);
for i=1:n
    I[i,i]=1;
end
Lf=L+I;
LL=inv(Lf);
sco=zeros(n);
for i=1:n
    sco[i]=LL[i,i];
end
selc3=zeros(k);
for i=1:k
    cho=argmax(sco);
    selc3[i]=cho;
    sco[cho]=-10000000;
#    cho=argmin(phi);
#    selc2[i*2-1]=cho;
#    phi[cho]=0;
end
L=lap(G);
for i=1:k
    L=lap(G);
    print(i,' ',lmd_LG(L,selc1[1:i],n));
    L=lap(G);
    println(' ',lmd_LG(L,selc2[1:i],n));
    #L=lap(G);
    #println(' ',lmd_LG(L,selc3[1:i],n));
end
close(fout)
close(fname)
#=
eigen(L, nev=1, which=:LR, v0 = ones(n))
B=zeros(n,k);
for i=1:k
    B[i,i]=1;
end
L=lap(G);
A=zeros(n,n);
for i=1:n
    for j=1:n
        A[i,j]=-L[i,j];
    end
end
for i=1:n
    A[i,i]=0;
end

D=zeros(n,n);
for i=1:n
    D[i,i]=L[i,i];
end

for i=1:4
    A[i*100,i*100]=1;
    D[i*100,i*100]+=1;
end

A=D^(-1/2)*A*D^(-1/2);
T=1000;
W1=B'*B;
W2=B*B';
for i=1:T
    W1=W1+B'*(A^(i))'*A^(i)*B;
    W1=(W1+W1')/2;
    f1=eigen(W1);
    W2=W2+A^(i)*B*B'*(A^(i))';
    W2=(W2+W2')/2;
    f2=eigen(W2);
    #=
    l=0;
    for j=1:n
        l+=WT[j,j];
    end
    =#
    #println(i," ",l);
    #=
    lmd=zeros(k);
    lmd=fi.values[n-k+1:n];
    for j=1:k
        lmd[j]=abs(lmd[j]);
    end
    =#
    println(i," ",f1.values[1]," ",f2.values[n-k+1]);
end
=#
#=
    L=lap(G);
    lead=argmax(L)[1];
    S=union(1:n);
    setdiff!(S,lead);
    Lr=L[S,S];
    f1=eigen(Lr);
    u=f1.vectors[:,1];
    if u[1]<0
        u=-u;
    end
    selc1=zeros(k);
    for i=1:k
        selc1[i]=argmax(u);
        u[argmax(u)]=0;
    end
    d=zeros(n-1);
    for i=1:n-1
        d[i]=Lr[i,i];
    end
    selc2=zeros(k);
    for i=1:k
        selc2[i]=argmax(d);
        d[argmax(d)]=0;
    end
    for i=1:k
        Lr=L[S,S];
        print(i,' ',lmdedge(Lr,selc1[1:i],i));
        Lr=L[S,S];
        println(' ',lmdedge(Lr,selc2[1:i],i));
    end
    =#
#=
L=lap(G);
n=G.n;
k=10;
nlmd=zeros(n);
L=lap(G);
D=zeros(n,n);
for i=1:n
    D[i,i]=L[i,i]^(-0.5);
end
L=D*L*D;
for i=1:n
    s=union(1:n);
    setdiff!(s,i);
    Li=L[s,s];
#    for j=1:n
#        Li[i,j]=0;
#        Li[j,i]=0;
#    end
    fi=eigen(Li);
    nlmd[i]=fi.values[1];
end

x=argmax(nlmd);
#L1 第一个点是精确的，后面选k-1个特征值最大的点
k=10;
selc1=zeros(k);
selc1[1]=x;
s=union(1:n);
setdiff!(s,x);
Li=L[s,s];
fi=eigen(Li);
#fi=lmdmin_val(Li,n-1);
u=fi.vectors[:,1];
if u[1]<0
    u=-u;
end
for i=2:k
    tmp=argmax(u);
    selc1[i]=tmp;
    if selc1[i]>=x
        selc1[i]+=1;
    end
    u[tmp]=0;
end


#=
dmax=maximum(L);
I=zeros(n,n);
for i=1:n
    I[i,i]=1;
end
A=dmax*I-L;
=#

#L2 第一个点是chenchen方法的，后面选k-1个特征值最大的点
L=lap(G);
A=L;
fa=eigen(A);
lmd=fa.values[n];
u=fa.vectors[:,n];
v=zeros(n);
for i=1:n
    v[i]=(2*lmd-A[i,i])*u[i]*u[i];
end
y=argmax(v);

k=10;
selc2=zeros(k);
selc2[1]=y;
s=union(1:n);
setdiff!(s,y);
L=lap(G);
Li=L[s,s];
fi=eigen(Li);
u=fi.vectors[:,1];
if u[1]<0
    u=-u;
end
for i=2:k
    tmp=argmax(u);
    selc2[i]=tmp;
    if selc2[i]>=y
        selc2[i]+=1;
    end
    u[tmp]=0;
end

#L3用chenchen方法选k个
L=lap(G);
dmax=maximum(L);
I=zeros(n,n);
for i=1:n
    I[i,i]=1;
end
A=dmax*I-L;

A=L;
fa=eigen(A);
lmd=fa.values[n];
u=fa.vectors[:,n];
v=zeros(n);
selc3=zeros(k);
for i=1:n
    v[i]=(2*lmd-A[i,i])*u[i]*u[i];
end
i=1;
selc3[1]=argmax(v);
for i=2:k
    S=union([]);
    for j=1:i-1
        union!(S,Int(selc3[j]));
    end
    B=A[:,S];
    b=B*u[S];
    scor=zeros(n);
    for j=1:n
        if j in S
            scor[j]=-1;
        else
            scor[j]=v[j]-2*b[j]*u[j];
        end
    end
    selc3[i]=argmax(scor);
end

#####4度最大k个
L=lap(G);
selc4=zeros(k);
deg=zeros(n);
for i=1:n
    deg[i]=L[i,i];
end
for i=1:k
       selc4[i]=argmax(deg);
       deg[argmax(deg)]=0;
   end
##########结果比较
L=lap(G);
for i=1:k
    println(i,' ',lmd_LG(L,selc1[1:i],n),' ',lmd_LG(L,selc2[1:i],n),' ',lmd_LG(L,selc3[1:i],n),' ',lmd_LG(L,selc4[1:i],n));
end
#=

LL=L;
for i=1:n
LL[x,i]=0;
LL[i,x]=0;
end
ff=eigen(LL);
vv=ff.vectors[:,2];
y=argmax(vv);
xx=0;yy=0;maxx=0;
corr=zeros(n,n);
for i=1:n-1
    for j=i+1:n
        s=union(1:n);
        setdiff!(s,i);
        setdiff!(s,j)
        Lij=L[s,s];
        fij=eigen(Lij);
        corr[i,j]=corr[j,i]=fij.values[1];
        #if corr[i,j]>maxx
        #    maxx=corr[i,j];
        #    xx=i;yy=j;
        #end
    end
end
=#
#=
    L=lap(G);
    n=G.n;
    lmin=zeros(n);
    fact=zeros(n);
    for i=1:n
        s=union(1:n);
        setdiff!(s,i);
        Li=L[s,s];
        fi=eigen(Li);
        fact[i]=fi.values[1];
        for j=1:n-1
            if sum(Li[j,:])>0
                Li[j,j]-=1;
            end
        end
        fi=eigen(Li);
        if abs(fi.values[2])<1e-5
            v=fi.values[3];
        else
            v=fi.values[2];
        end
        lmin[i]=L[i,i]/(n-1)*(1-2*sqrt(L[i,i])/v);
    end
    for i=1:n
        println(i,' ',lmin[i],' ',fact[i],' ',L[i,i]);
    end

L=lap(G);
f=eigen(L);
u2=f.vectors[:,2];
B=zeros(n,n);
for i=1:n
    for j=1:n
        B[i,j]=-L[i,j]*u2[i]*u2[j];
    end
end
for i=1:n
    B[i,i]=0;
end
infl=zeros(n);
for i=1:n
    infl[i]=sum(B[i,:])-B[i,i];
end
=#
#=
L=lap(G);
n=G.n;
f=eigen(L);
lambda=f.values[2];
u0=f.vectors[:,2];
incr=zeros(n);
fact=zeros(n);
for i=1:n
    s=union(1:n);
    setdiff!(s,i);
    LL=L[s,s];
    fi=eigen(LL);
    u=fi.vectors[:,1];
    if u[1]<0
        u=-u;
    end
    if u[2]<0
        u=-u;
    end
    if u[3]<0
        u=-u;
    end
    II=zeros(n-1,n-1);
    for j=1:n-1
        II[j,j]=1;
    end
    #LL=(eigen(LL).values[n-1]+1)*II-LL;
    #LL=II-0.5*LL;
    v=ones(n-1);
    #v[1:i-1]=u0[1:i-1];
    #v[i:n-1]=u0[i+1:n];

    for j=1:n-1
        #v[j]=LL[j,j];
    end
    LL=(maximum(LL)+1)*II-LL;
    #LL=II-LL;
    #v=v./norm(v);
    v=LL*v;
    v=v./norm(v);
    E=spzeros(n,n);
    for j=1:n
        E[i,j]=-L[i,j];
        E[j,i]=-L[j,i];
    end
    #E=-E;
    u1=u0;
    u2=zeros(n);
    u2[1:i-1]=v[1:i-1];
    u2[i+1:n]=v[i:n-1];
u1=ones(n);
incr[i]=(u1'*E*u2)/(u1'*u2);

#    incr[i]=(u1'*E*u2)/(u1'*u2);



    fact[i]=fi.values[1];
    u3=zeros(n);
    u3[1:i-1]=fi.vectors[1:i-1,1];
    u3[i+1:n]=fi.vectors[i:n-1,1];
    if u3[1]<0
        u3=-u3;
    end
    if u3[2]<0
        u3=-u3;
    end
    if u3[3]<0
        u3=-u3;
    end
    incc=(u1'*E*u3)/(u1'*u3);
    #println(i,' ',norm(u2-u3));
    #println(i,' ',incr[i],' ',fact[i]-lambda,' ',L[i,i]);
end

x=argmax(incr);
y=argmin(incr);
z=argmax(fact);
println(x,' ',y,' ',z);
println(fact[z]/fact[x],' ',fact[z]/fact[y]);
##############
=#
#=
I=zeros(n,n);
for i=1:n
   I[i,i]=1;
end
delta=0.01;
C=I-delta*L-1/n*ones(n,n);
fc=eigen(C);
lc=fc.values[n];
lmd=(1-lc)/delta;
s=1-delta*lmd+delta;
p=delta*lmd*lmd-lmd;
q=1-delta*lmd;
uu=fc.vectors[:,n];
u=eigen(L).vectors[:,2];
old=zeros(n);
for i=1:n
    E=spzeros(n,n);
    for j=1:n
        E[i,j]=-L[i,j];
        E[j,i]=-L[j,i];
    end
    t=u'*E*u;
    old[i]=(s*t+p*u[i])/(q+delta*t-u[i]^2);
end
for i=1:n
    #println(i,' ',old[i],' ',fact[i]-lambda);
end
x=argmax(incr);
y=argmin(incr);
z=argmax(fact);
println(x,' ',y,' ',z);
println(fact[z]/fact[x],' ',fact[z]/fact[y]);
println(123);

=#

#    parB=getB(G,G.m);
#    parL=lapsp(G);
#for i=1:G.n
#        parL[i,i]=0;
#    end
#    parL=-parL;
#  f = approxCholLap(parL, tol=1e-8);
#      m = size(parB, 1);
#      n = size(parL, 1);
#      k = 10*round(Int, log(n)); # number of dims for JL
#      rst = 0;
#      for i = 1:k
#        r = rand([1,-1],m);
#        v = r'*parB;
#        y=f(v');
#        rst+=y'*y;
#      end
#println(rst/k);
#L=lap(G);
#invL=inv(L-1/G.n*ones(G.n,G.n))+1/G.n*ones(G.n,G.n);
#tr=0;
#for i=1:G.n
#    tr+=invL[i,i];
#end
#println(tr);
#close(fname)
    #fout = open(string("data/","finnal.ans"),"a")
    #println("Now running file:",str[1],"with eps=",ep);
#    for k= 1 : 3
#        ep = 0.1*k;
#        println("Now running file:",str[1],"with eps=",ep);
#        fout = open(string("data/","time.ans"),"a");
#        print(fout,str[1]," ");
#        greed = gre(G,G.k ,ep,eta,fout);
#        println("gre=",greed);
#        println(fout," ");
#        close(fout);
#    end
#    fout = open(string("data/","time.ans"),"a");
#    print(fout,str[1]," ");
#    println("Now running file:",str[1]," exactly");
#    exact=exa(G,G.k,eta,fout);
#    println("exa=",exact);
#    close(fout);
#end

#close(fname)
=#
