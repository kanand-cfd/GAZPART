include Makefile.imft


MODOBJS	= \
mod_dns_dim.o \
mod_check_cpu.o \
mod_cputime.o \
mod_ensight_var.o \
mod_fluid_variable.o \
mod_forcing.o \
mod_geometric_variable.o \
mod_param_phys.o\
mod_quaternion.o \
mod_particle_parallel.o \
mod_scalar_variable.o \
mod_statistics.o \
mod_work_arrays.o \
mod_collision.o \
mod_random_var.o \
MPI_structures.o \
utilities_cpu.o \
exchange_part.o \
minpack_reduc.o


FLUOBJS	= \
add_forcing3.o \
adv_forcing.o \
adv_fluid.o \
allocate_arrays.o\
nonlinear_terms.o \
check_cfl.o \
cputime_info.o \
corrtwopts.o \
divergence.o \
fluid_gradient.o \
inforun.o\
init_run.o \
initiation_fluid.o\
initiation_forcing_new.o\
meshing.o \
openclose.o \
proj_divfree.o \
print_timestat.o \
print_lagfunction.o \
readparam.o \
save_fluid.o\
save_forcing.o\
save_forcing_new.o\
spec3d.o\
stat_fluid.o \
tremain.o \
main.o


PAROBJS	= \
adv_particle.o \
boundary_particle.o \
boundary_particle_single.o \
boundary_particle_yz.o\
collision.o \
collision_box_connectivity.o \
collision_detection_connectivity.o\
collision_initiation.o \
collision_interaction.o\
initiation_particle_position.o \
initiation_particle_velocity.o \
initiation_frozen_flow.o \
interpH.o \
interp_lag1.o\
interp_lag2.o\
interp_lag3.o\
interp_lag4.o \
frozen_interp_lag3.o \
wall_interp_lag3.o\
print_particle.o \
stat_particle.o \
save_particle.o \
time_averaged_stat.o \
stat_part_collision.o \
init_channel_stat.o \
locate_part.o\
stat_part_channel.o\
print_timestat_channel.o\
before_coll.o\
after_coll.o\
stat_coll_channel.o

ELLOBJS = \
adv_particle_new.o \
adv_quaternion.o \
adv_angular_velocity.o \
adv_angular_velocity_fluid.o \
adv_quat_angvel_fluid.o \
invert_ndim3_matrix.o \
calculate_ellipsoid_matrix.o \
ellipsoid_pair_depth.o \
ellipsoid_contact_detection.o \
collision_response.o \
transform_basis.o \
fluid_force_ellipsoid.o \
fluid_torque_ellipsoid.o \
wall_boundary_condition.o \
closest_point_plane.o \
max_encroachment_ellipsoid.o \
wall_ellipsoid_rebound.o \
stat_partfluid.o \
print_colliding_pair.o \
ellipsoid_wall_boundary_new.o 


MOD	= *.mod
OBJS	= $(MODOBJS) $(FLUOBJS) $(PAROBJS) $(ELLOBJS)

EXEC    =  ../../ell_gaspart 


$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) $(LIBS) $(INCLUDE) -o $@

clean:
	$(RM) $(OBJS) $(MOD) $(EXEC)

cpart:
	$(RM) $(PAROBJS) $(MOD) $(EXEC)

.SUFFIXES: .o .f90 

.f90.o:
	$(FC)  -c $(FFLAGS) $(INCLUDE) $<
