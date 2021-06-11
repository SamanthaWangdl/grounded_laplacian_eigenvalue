using LinearAlgebra
using SparseArrays
using Laplacians
using Random

function lmdmin_val(L,n)
	I=spzeros(n,n);
	for i=1:n
		I[i,i]=n;
	end
	for i=1:n
		for j=1:n
			L[i,j]=I[i,j]-L[i,j];
		end
	end
	u=randn(n);
	u=u./norm(u);
	y=L*u;
	y=y./norm(y);
	while norm(y-u)>1e-2
		u=y;
		y=L*u;
		y=y./norm(y);
	end
	return n-u'*L*y;
end

function lmdedge(L,s,k)
	for i=1:k
		L[Int(s[i]),Int(s[i])]+=1;
	end
	f=eigen(L);
	return f.values[1];
end

function lmd_LG(L,s,n)
	S=union(1:n);
	setdiff!(S,s);
	LG=L[S,S];
	f=eigvals(LG);
	k=1;
	while f[k]<1e-8
		k+=1;
	end
	return f[k];
end

function lap(G :: Graph)
    F = zeros(G.n, G.n);
    for i = 1 : G.m
        F[G.u[i], G.v[i]] -= 1
        F[G.v[i], G.u[i]] -= 1
        F[G.u[i], G.u[i]] += 1
        F[G.v[i], G.v[i]] += 1
    end
    return F
end

function lapsp(G :: Graph)
    F = spzeros(G.n, G.n)
    for i = 1 : G.m
        F[G.u[i], G.v[i]] -= 1
        F[G.v[i], G.u[i]] -= 1
        F[G.u[i], G.u[i]] += 1
        F[G.v[i], G.v[i]] += 1
    end
    return F
end

function getB(G :: Graph, mm)
    B :: SparseMatrixCSC{Float64} = spzeros(mm,G.n)
	jishu=0;
    for i = 1 : G.m
    #    if (G.u[i]>G.n0+G.n1 && G.v[i]>G.n0+G.n1)
			jishu+=1;
            B[jishu, G.u[i]] = 1
            B[jishu, G.v[i]] = -1
    #    end
    end
    return B
end

function getD(G :: Graph)
    D = spzeros(G.n-G.n0-G.n1,G.n-G.n0-G.n1);
    for i = 1 : G.m
        if (G.u[i]<=G.n0+G.n1) && (G.v[i]>G.n0+G.n1)
            D[G.v[i]-G.n0-G.n1] += 1;
        end
        if (G.v[i]<=G.n0+G.n1) && (G.u[i]>G.n0+G.n1)
            D[G.u[i]-G.n0-G.n1] += 1;
        end
    end
    return D;
end

function calc(G :: Graph ,i ,Linv, sel)
    e=spzeros(G.n-G.n0-G.n1,1)
    e[i,1]=1
    vect = Linv*e-Linv*e*(e'*Linv*sel);
    ans = sum(vect)/(1+Linv[i,i]);
    return ans
end


function findmaxexa(G :: Graph, Linv ,sel,eta)
    bestv = 0;
    tmp = 0
    maxinc = 0
    for i = G.n0+G.n1+1 : G.n
        if sel[i-(G.n0+G.n1)]< G.n1
            tmp = calc(G, i-G.n0-G.n1 , Linv ,sel)
            if tmp > maxinc && (Linv*sel)[i-G.n0-G.n1]<eta
                maxinc = tmp
                bestv = i
            end
			#print(tmp," ")
        end
    end
    return bestv
end


function exa(G :: Graph,k,eta,fout)
    L = lap(G);
    sel = zeros(G.n - G.n0- G.n1);
    for i = G.n0+G.n1+1 : G.n
        sel[i-(G.n0+G.n1)] =  sum(L[i,1:G.n0]) - sum(L[i,1:G.n0+G.n1]);
    end
	T1=time();
    Linv = inv(L[G.n0+G.n1+1:G.n,G.n0+G.n1+1:G.n]);
	firans = sum(Linv * sel);
	j=0;
	T2=time();
	#println(" ");
    for i = 1 : k
        j = findmaxexa(G, Linv ,sel,eta);
		#print(j-G.n0-G.n1," ")
        sel[j-G.n0-G.n1] += 1;
        e=spzeros(G.n-G.n0-G.n1,1);
        e[j-G.n0-G.n1,1]=1;
        Linv = Linv - Linv*e*e'*Linv/(1+Linv[j-G.n0-G.n1,j-G.n0-G.n1]);
    end
	#println("");
	T3=time();
	println(fout,"exact basetime=",T2-T1,",",k," times=",T3-T2);
    ansv = Linv * sel;
    return sum(ansv);
end

function gre(G :: Graph, k ,ep,eta,fout)
	finans= 0 ;
    d0=ep/(3*sqrt(G.n));
	#d0=eps/(3);
    #d=4*eps/(3*n^2);
	d=4*ep/(3*G.n^2)/(G.n^2);
	JLfac=1.0;
    #t= round(Int64,12*12*24*log2(G.n)/eps^2);
	t= round(Int64,JLfac*log(G.n)/(ep^2))+1;
	L = lapsp(G);
	nn = G.n-G.n0-G.n1;
    LL = L[G.n0+G.n1+1:G.n,G.n0+G.n1+1:G.n];
    mm = 0;
    for i = 1: G.m
        if ((G.u[i]>G.n0+G.n1) && (G.v[i]>G.n0+G.n1))
            mm += 1;
        end
    end
	B = getB(G,mm);
    D = getD(G);
	Dadd = spzeros(G.n-G.n0-G.n1,G.n-G.n0-G.n1);
    #sqrt.(D);
    rng = MersenneTwister(Int(round(time() * 10000.0)));
    #f = approxchol_sddm(LL, tol=d0)

    sel = zeros(nn);
    for i = G.n0+G.n1+1 : G.n
        sel[i-(G.n0+G.n1)] =  sum(L[i,1:G.n0]) - sum(L[i,1:G.n0+G.n1]);
    end
    o1 = ones(nn);
	maxx = 0;
	z1 = zeros(nn, 1)
	z2 = zeros(nn, 1)
	f = approxchol_sddm(LL, tol=1e-8);
	firans =sum(f(sel));
	#r2 = randn(rng, nn, 1);
	#r1 = randn(rng, mm, 1);

	T1=time();

    for i = 1 : k

##################


####f = Laplacians.approxchol_sddm(LL+Dadd, tol=d0);
f = Laplacians.approxchol_sddm(LL, tol=d0,verbose=false);
        h1=f(o1);
        h2=f(sel);

#这里就是用求解器的步骤，f = approxchol_sddm(L, tol=d0);L是拉普拉斯矩阵，tol=d0是误差，选取1e-6完全够了，
#h1=f(o1)就是返回L的伪逆*o1向量
#####################



		Xe = zeros(nn);
		Ye = zeros(nn);
        for j = 1 : t
            r2 = randn(rng, nn, 1);
			r1 = randn(rng, mm, 1);
			xx = B'*r1;
			yy = sqrt.(D)*r2;
            f = approxchol_sddm(LL, tol=d);
			z1[:, 1] = f(xx[:, 1]);
			z2[:, 1] = f(yy[:, 1]);
			for p = 1 : nn
				Xe[p] += z1[p,1]^2;
				Ye[p] += z2[p,1]^2;
			end
		end
		add = zeros(G.n-G.n0-G.n1);
		for j = 1 : G.n-G.n0-G.n1
			add[j] += h1[j]*(1-h2[j])/(1+Xe[j]/t+Ye[j]/t);
		end
		maxx=0;
		bestv=0;
		for j = 1 : nn
			if (sel[j]<G.n1 && add[j]>maxx && h2[j]<eta)
				maxx = add[j]
				bestv = j
			end
		end
		#println(" ");
		#for j =1 : G.n-G.n0-G.n1
		#	print(add[j]," ")
		#end
		sel[bestv]+=1;
		#LL[bestv,bestv]+=1;
		D[bestv,bestv]+=1;
		Dadd[bestv,bestv]+=1;
		finans += maxx;
		#print(bestv," ")
	end
	T2=time();
	#println(" ");
	f = approxchol_sddm(LL+Dadd, tol=1e-8);
	finans =sum(f(sel));
	println(fout,"greed ",k," times=",T2-T1," eps=",ep);
	return finans;
end



function appxInvTrace(parL::SparseMatrixCSC{Float64}, parB::SparseMatrixCSC{Float64})

  f = approxchol_sddm(parL,tol=1e-6, verbose=false)
  m = size(parB, 1)
  n = size(parL, 1)
  k = round(Int, log(n)) # number of dims for JL
  print(n)
  rst = 0

  for i = 1:k
  #  V[:,i] .= f(UR[:,i])
    r = rand([1,-1],m,1)
    #ur = U'*r
    v = r'*parB
    print(v)
    v = f(v')
    #rst += sum(sqrt(v.*v))
  end

  return rst/k;

end
