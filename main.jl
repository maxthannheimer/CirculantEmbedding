include("functions.jl")





#define grid, simulate and plot

gridsize=1000 #size of grid is gridsize*gridsize, simulated grid is (2grid_size)²  #
param=[1.0,1.0] # param[1]=c, param[2]=alpha
#param defines the cavariance function: 
#cov(x,y)=      c ⋅ (     ||x||^α + ||y||^α − ||x-y||^α       )
num_sim=1 #number of simulations

#FBM Simulation:
@time field=FBM_simu_fast(param, gridsize, num_sim)[1]
#result is matrix representing gridpoints points in [0,1] × [0,1] 
#in the upper left corner is (0,0)
# x_axis grows to the right (columnwise), upper right corner is (1,0)
# y_axis grows downswoards (rowwise), lower left corner is (0,1)


#grid points in x and y 
tx=ty=(1:gridsize)/gridsize

plotd=surface(tx,ty,field,title="$grid_size × $grid_size FBM simulation, α=$(param[2]), c= $(param[1])")
