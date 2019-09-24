#import pdp_lr_b as d
import laparcgis7 as d
#import laparcgis6geo_tested as d
import time,sys
seed=-1
psize=10
t=100
spp=0
n=4
#fn="C:\\PDP\\fldp_second\\gy2s_20c.txt"
#fn2="C:\\PDP\\fldp_second\\gy2s_connectivity.txt"
#fn="C:\\PDP\\fldp_second\\gys_37c.txt"
#fn2="C:\\PDP\\fldp_second\\gys_connectivity.txt"
fn="C:\\PDP\\fldp_second\\arcgis\\zys4_units.txt"
fn2="C:\\PDP\\fldp_second\\arcgis\\zys_connectivity.txt"

t0=time.time()
if len(sys.argv) >= 2:
    fn=sys.argv[1]
if len(sys.argv) >= 3:
    fn2=sys.argv[2]

if len(sys.argv) >= 4:
    n=int(sys.argv[3])
if len(sys.argv) >= 5:
    psize=int(sys.argv[4])
if len(sys.argv) >= 6:
    t=int(sys.argv[5])
if len(sys.argv) >= 7:
    spp=int(sys.argv[6])
if len(sys.argv) >= 8:
    seed=int(sys.argv[7])

d.solver_message=1
d.readfile(fn,fn2)
#d.read_bm_instance(fn)

d.location_problem=2   #0 ap, sap; 1 flp; 2 pdp
d.pop_dis_coeff=1.0
d.pop_deviation=0.0

d.fixed_cost_obj=0
d.spatial_contiguity=01
d.allowing_less_facilities=0
d.initial_solution_method=0 #0greedy,1 lp, 2 ap, 9 LR
#for PDP, 0=KM, 1,
d.mip_solver="cplex"
d.operators_selected=[0,1]
d.solution_similarity_limit=3#max(10.0,100.0/n)
d.solver_message=01

d.acceptanceRule="ils"
d.heuristic_time_limit=t
d.spp_modeling=spp
d.multi_start_count=psize
d.seed=seed
d.max_num_facility=n

d.initialize_instance()

#d.solve()
#d.max_num_facility=n

#d.initialize_instance()
#d.initial_solution(0)
#print d.centersID
#d.mipmodel_pulp() 

#d.ils_lr(n,psize,t,spp,seed)
d.ils(n,psize,t,spp,seed)
###d.ils_lr_pdp(n,psize,t,spp,seed)
#d.mip(n,2,1,t)
#d.print_solution()

print "final best:",d.biobjective,d.objective,d.objective_overload
print time.time()-t0
print len(set(d.centersID)),d.centersID
for x in d.district_info: 
    if x[0]>0: print x
print d.avg_dis_min