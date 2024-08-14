#import Pkg; Pkg.add("Distributions")
using LinearAlgebra
using Random, Distributions
using Plots
using FFTW


#################################################
#################################################
#Simulation of fractinal Brownian Motion in 2d via Circulant Embedding with complexity O(n log(n) ) for n gridpoints
#From C.R. Dietrich & G.N. Newsam (1997):
#        Fast and exact simulation of stationary Gaussian processes through circulant embedding of the covariance matrix.

#via algorithm from Volker Schmidt's Book:
#Stochastic Geometry, Spatial Statistics and Random Fields -- Models and Algorithms
# Chapter: 12.4.4 Fractional Brownian Field



#################################################
#################################################

#create functions for embedding simulation, one for 2dfft, one for modified cov function and one for the simulation

#mathlab fft2 analog function
function fft2(A)
    # Apply 1D FFT along each dimension
    fft_rows = fft(A, 1)
    fft_cols = fft(fft_rows, 2)
    # Return the result
    return fft_cols
end

#by hand calculation of fourier transform
function getFourierMatrix(n)
    F=rand(n,n)*im
    for i in 0:(n-1),j in 0:(n-1)
        F[i+1,j+1]=exp(2*pi*im*i*j/n)
    end
    F   
end

function fft2_hardcode(A)
    F=getFourierMatrix(size(A,1))
    fft_rows = F*A
    fft_cols = (F*fft_rows')'
    # Return the result
    return fft_cols
end


#modified cov function
function rho(x,y,R,alpha)
    #embedding of covariance function on a larger [0,R] × [0,R] Grid
    if alpha<=1.1 #alpha=2 Hurst param
        beta=0
        c_2=alpha/2
        c_0=1-alpha/2
    else
        beta=alpha*(2-alpha)/(3*R*(R^2-1))
        c_2=(alpha-beta*(R-1)^2*(R+2))/2
        c_0=beta*(R-1)^3+1-c_2
    end
    #create cont isotropic cov function
    r=sqrt((x[1]-y[1])^2+(x[2]-y[2])^2)
    if r<=1
        out=c_0-r^alpha+c_2*r^2
    elseif R==1 #caution?
        out=0
    elseif r<=R
        out=beta*(R-r)^3/r
    else
        out=0
    end
    return (out,c_0,c_2)
end

function FBM_simu_fast(param,gridsize,num_sim) 
    alpha=param[2]
    c_sqrt=sqrt(param[1])
    H=alpha/2 #Hurst Param
    if alpha<=1.1
        R=1
    else
        R=2 #expanded gridsize, region of interest is [0,1]×[0.1]
    end
    n=m=R*gridsize #size of grid is m ∗ n, cov mat size is n^2*m^2
    tx=(1:n)/n*R; ty=(1:m)/m*R #grid points in x and y 
    Rows=zeros(m,n);
    for i in 1:n 
        for j in 1:m
            Rows[j,i]=rho([tx[i],ty[j]],[tx[1],ty[1]],R,2*H)[1]
        end
    end
    BlkCirc_row=[Rows  Rows[:,end-1:-1:2] ;  Rows[end-1:-1:2,:]  Rows[end-1:-1:2,end-1:-1:2 ]  ]
    #calculate eigenvalues via fft
    eig_vals=real(fft2(BlkCirc_row)/(4*(m-1)*(n-1)))
    #optional:
    #set small values to zero:
    #eps = 10^-8
    #eig_vals=eig_vals.*(abs.(eig_vals) .> eps)
    #eig_vals=eig_vals.*(eig_vals .>= 0)
    eig_vals=sqrt.(eig_vals)
   
    res1=[Matrix{Float64}(undef,gridsize,gridsize) for i in 1:num_sim]
    #one can get two times as many simulations for free by using the imaginary and real part 
    #of the complex gaussian, but they are dependend
    
    #res2=[Matrix{Float64}(undef,gridsize,gridsize) for i in 1:num_sim]
    
    for trial in 1:num_sim
        #generate field with covariance given by block circulant matrix
        Z= randn(2*(m-1),2*(n-1)) + im* randn(2*(m-1),2*(n-1)) 
        #fft to calculate diag*Z
        F=fft2(eig_vals.*Z)
        #extract subblock with desired cov variance
        F=F[1:gridsize,1:gridsize]
        (out,c_0,c_2)=rho([0,0],[0,0],R,2*H)
        #optional two (dependend) real fields
        field1=real(F)
        #field2=imag(F)
        #set field zero at origin
        field1=field1.-field1[1,1]
        #field2=field2.-field2[1,1]
        #make correction for embedding with a term c_2*r^2
        X_grid_comp = [i for i in tx[1:gridsize], j in tx[1:gridsize]].-tx[1]
        Y_grid_comp = [j for i in tx[1:gridsize], j in tx[1:gridsize]].-tx[1]
    
        res1[trial]=field1 + (X_grid_comp*randn()+Y_grid_comp*randn() ) .*sqrt(2*c_2)
        #res2[trial]=field2 +  (X_grid_comp*randn()+Y_grid_comp*randn() ) .*sqrt(2*c_2)
    end

    #(c_sqrt.*res1,c_sqrt.*res2)
    c_sqrt.*res1
end

