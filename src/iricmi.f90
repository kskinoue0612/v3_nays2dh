module iricmi
	implicit none

	interface iricmi_read_grid2d_integer_node
		module procedure iricmi_read_grid2d_integer_node_1d
		module procedure iricmi_read_grid2d_integer_node_2d
	end interface

	interface iricmi_read_grid2d_real_node
		module procedure iricmi_read_grid2d_real_node_1d
		module procedure iricmi_read_grid2d_real_node_2d
	end interface

	interface iricmi_read_grid2d_complex_node
		module procedure iricmi_read_grid2d_complex_node_1d
		module procedure iricmi_read_grid2d_complex_node_2d
	end interface

	interface iricmi_read_grid2d_integer_cell
		module procedure iricmi_read_grid2d_integer_cell_1d
		module procedure iricmi_read_grid2d_integer_cell_2d
	end interface

	interface iricmi_read_grid2d_real_cell
		module procedure iricmi_read_grid2d_real_cell_1d
		module procedure iricmi_read_grid2d_real_cell_2d
	end interface

	interface iricmi_read_grid2d_complex_cell
		module procedure iricmi_read_grid2d_complex_cell_1d
		module procedure iricmi_read_grid2d_complex_cell_2d
	end interface

	interface iricmi_rout_grid2d_integer_node
		module procedure iricmi_rout_grid2d_integer_node_1d
		module procedure iricmi_rout_grid2d_integer_node_2d
	end interface

	interface iricmi_rout_grid2d_real_node
		module procedure iricmi_rout_grid2d_real_node_1d
		module procedure iricmi_rout_grid2d_real_node_2d
	end interface

	interface iricmi_rout_grid2d_integer_cell
		module procedure iricmi_rout_grid2d_integer_cell_1d
		module procedure iricmi_rout_grid2d_integer_cell_2d
	end interface

	interface iricmi_rout_grid2d_real_cell
		module procedure iricmi_rout_grid2d_real_cell_1d
		module procedure iricmi_rout_grid2d_real_cell_2d
	end interface

contains

	subroutine iricmi_check_cancel(canceled, ier)
		integer, intent(out):: canceled
		integer, intent(out):: ier

		call iricmi_check_cancel_f2c(canceled, ier)

	end subroutine

	subroutine iricmi_model_init(ier)
		integer, intent(out):: ier

		call iricmi_model_init_f2c(ier)

	end subroutine

	subroutine iricmi_model_terminate(ier)
		integer, intent(out):: ier

		call iricmi_model_terminate_f2c(ier)

	end subroutine

	subroutine iricmi_model_sync(ier)
		integer, intent(out):: ier

		call iricmi_model_sync_f2c(ier)

	end subroutine

	subroutine iricmi_model_dump(ier)
		integer, intent(out):: ier

		call iricmi_model_dump_f2c(ier)

	end subroutine

	subroutine iricmi_calc_init(ier)
		integer, intent(out):: ier

		call iricmi_calc_init_f2c(ier)

	end subroutine

	subroutine iricmi_calc_terminate(ier)
		integer, intent(out):: ier

		call iricmi_calc_terminate_f2c(ier)

	end subroutine

	subroutine iricmi_calc_sync_receive(ier)
		integer, intent(out):: ier

		call iricmi_calc_sync_receive_f2c(ier)

	end subroutine

	subroutine iricmi_calc_sync_send(ier)
		integer, intent(out):: ier

		call iricmi_calc_sync_send_f2c(ier)

	end subroutine

	subroutine iricmi_calc_dump(ier)
		integer, intent(out):: ier

		call iricmi_calc_dump_f2c(ier)

	end subroutine

	subroutine iricmi_error_print()

		call iricmi_error_print_f2c()

	end subroutine

	subroutine iricmi_error_print_and_exit()

		call iricmi_error_print_and_exit_f2c()

	end subroutine

	subroutine iricmi_read_integer(name, val, ier)
		character(*), intent(in):: name
		integer, intent(out):: val
		integer, intent(out):: ier

		call iricmi_read_integer_f2c(name, val, ier)

	end subroutine

	subroutine iricmi_read_real(name, val, ier)
		character(*), intent(in):: name
		double precision, intent(out):: val
		integer, intent(out):: ier

		call iricmi_read_real_f2c(name, val, ier)

	end subroutine

	subroutine iricmi_read_realsingle(name, val, ier)
		character(*), intent(in):: name
		real, intent(out):: val
		integer, intent(out):: ier

		call iricmi_read_realsingle_f2c(name, val, ier)

	end subroutine

	subroutine iricmi_read_stringlen(name, len, ier)
		character(*), intent(in):: name
		integer, intent(out):: len
		integer, intent(out):: ier

		call iricmi_read_stringlen_f2c(name, len, ier)

	end subroutine

	subroutine iricmi_read_string(name, str, ier)
		character(*), intent(in):: name
		character(*), intent(out):: str
		integer, intent(out):: ier

		call iricmi_read_string_f2c(name, str, ier)

	end subroutine

	subroutine iricmi_read_functional_size(name, size, ier)
		character(*), intent(in):: name
		integer, intent(out):: size
		integer, intent(out):: ier

		call iricmi_read_functional_size_f2c(name, size, ier)

	end subroutine

	subroutine iricmi_read_functional_vals(name, x, y, ier)
		character(*), intent(in):: name
		double precision, dimension(:), allocatable, intent(out):: x
		double precision, dimension(:), allocatable, intent(out):: y
		integer, intent(out):: ier

		call iricmi_read_functional_vals_f2c(name, x, y, ier)

	end subroutine

	subroutine iricmi_read_functional_valwithname(name, paramname, val, ier)
		character(*), intent(in):: name
		character(*), intent(in):: paramname
		double precision, dimension(:), allocatable, intent(out):: val
		integer, intent(out):: ier

		call iricmi_read_functional_valwithname_f2c(name, paramname, val, ier)

	end subroutine

	subroutine iricmi_read_functional_vals_realsingle(name, x, y, ier)
		character(*), intent(in):: name
		real, dimension(:), allocatable, intent(out):: x
		real, dimension(:), allocatable, intent(out):: y
		integer, intent(out):: ier

		call iricmi_read_functional_vals_realsingle_f2c(name, x, y, ier)

	end subroutine

	subroutine iricmi_read_functional_valwithname_realsingle(name, paramname, val, ier)
		character(*), intent(in):: name
		character(*), intent(in):: paramname
		real, dimension(:), allocatable, intent(out):: val
		integer, intent(out):: ier

		call iricmi_read_functional_valwithname_realsingle_f2c(name, paramname, val, ier)

	end subroutine

	subroutine iricmi_read_grid_node_count(count, ier)
		integer, intent(out):: count
		integer, intent(out):: ier

		call iricmi_read_grid_node_count_f2c(count, ier)

	end subroutine

	subroutine iricmi_read_grid_cell_count(count, ier)
		integer, intent(out):: count
		integer, intent(out):: ier

		call iricmi_read_grid_cell_count_f2c(count, ier)

	end subroutine

	subroutine iricmi_read_grid2d_str_size(isize, jsize, ier)
		integer, intent(out):: isize
		integer, intent(out):: jsize
		integer, intent(out):: ier

		call iricmi_read_grid2d_str_size_f2c(isize, jsize, ier)

	end subroutine

	subroutine iricmi_read_grid2d_unstr_size(nodessize, ier)
		integer, intent(out):: nodessize
		integer, intent(out):: ier

		call iricmi_read_grid2d_unstr_size_f2c(nodessize, ier)

	end subroutine

	subroutine iricmi_read_grid2d_unstr_cellsize(cellssize, ier)
		integer, intent(out):: cellssize
		integer, intent(out):: ier

		call iricmi_read_grid2d_unstr_cellsize_f2c(cellssize, ier)

	end subroutine

	subroutine iricmi_read_grid2d_unstr_cells(cells, ier)
		integer, dimension(:), allocatable, intent(out):: cells
		integer, intent(out):: ier

		call iricmi_read_grid2d_unstr_cells_f2c(cells, ier)

	end subroutine

	subroutine iricmi_read_grid2d_coords(x, y, ier)
		double precision, dimension(:,:), intent(out):: x
		double precision, dimension(:,:), intent(out):: y
		integer, intent(out):: ier

		call iricmi_read_grid2d_coords_f2c(x, y, ier)

	end subroutine

	subroutine iricmi_read_grid2d_integer_node_1d(name, vals, ier)
		character(*), intent(in):: name
		integer, dimension(:), intent(out):: vals
		integer, intent(out):: ier

		call iricmi_read_grid2d_integer_node_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_read_grid2d_integer_node_2d(name, vals, ier)
		character(*), intent(in):: name
		integer, dimension(:,:), intent(out):: vals
		integer, intent(out):: ier

		call iricmi_read_grid2d_integer_node_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_read_grid2d_real_node_1d(name, vals, ier)
		character(*), intent(in):: name
		double precision, dimension(:), intent(out):: vals
		integer, intent(out):: ier

		call iricmi_read_grid2d_real_node_f2c(name, vals, ier)

	end subroutine
	subroutine iricmi_read_grid2d_real_node_2d(name, vals, ier)
		character(*), intent(in):: name
		double precision, dimension(:,:), intent(out):: vals
		integer, intent(out):: ier

		call iricmi_read_grid2d_real_node_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_read_grid2d_complex_node_1d(name, vals, ier)
		character(*), intent(in):: name
		integer, dimension(:), intent(out):: vals
		integer, intent(out):: ier

		call iricmi_read_grid2d_complex_node_f2c(name, vals, ier)

	end subroutine
	subroutine iricmi_read_grid2d_complex_node_2d(name, vals, ier)
		character(*), intent(in):: name
		integer, dimension(:,:), intent(out):: vals
		integer, intent(out):: ier

		call iricmi_read_grid2d_complex_node_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_read_grid2d_integer_cell_1d(name, vals, ier)
		character(*), intent(in):: name
		integer, dimension(:), intent(out):: vals
		integer, intent(out):: ier

		call iricmi_read_grid2d_integer_cell_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_read_grid2d_integer_cell_2d(name, vals, ier)
		character(*), intent(in):: name
		integer, dimension(:,:), intent(out):: vals
		integer, intent(out):: ier

		call iricmi_read_grid2d_integer_cell_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_read_grid2d_real_cell_1d(name, vals, ier)
		character(*), intent(in):: name
		double precision, dimension(:), intent(out):: vals
		integer, intent(out):: ier

		call iricmi_read_grid2d_real_cell_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_read_grid2d_real_cell_2d(name, vals, ier)
		character(*), intent(in):: name
		double precision, dimension(:,:), intent(out):: vals
		integer, intent(out):: ier

		call iricmi_read_grid2d_real_cell_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_read_grid2d_complex_cell_1d(name, vals, ier)
		character(*), intent(in):: name
		integer, dimension(:), intent(out):: vals
		integer, intent(out):: ier

		call iricmi_read_grid2d_complex_cell_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_read_grid2d_complex_cell_2d(name, vals, ier)
		character(*), intent(in):: name
		integer, dimension(:,:), intent(out):: vals
		integer, intent(out):: ier

		call iricmi_read_grid2d_complex_cell_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_read_bc_count(type, num, ier)
		character(*), intent(in):: type
		integer, intent(out):: num
		integer, intent(out):: ier

		call iricmi_read_bc_count_f2c(type, num, ier)

	end subroutine

	subroutine iricmi_read_bc_indicessize(type, num, size, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		integer, intent(out):: size
		integer, intent(out):: ier

		call iricmi_read_bc_indicessize_f2c(type, num, size, ier)

	end subroutine

	subroutine iricmi_read_bc_indices(type, num, indices, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		integer, dimension(:), intent(out):: indices
		integer, intent(out):: ier

		call iricmi_read_bc_indices_f2c(type, num, indices, ier)

	end subroutine

	subroutine iricmi_read_bc_integer(type, num, name, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		integer, intent(out):: val
		integer, intent(out):: ier

		call iricmi_read_bc_integer_f2c(type, num, name, val, ier)

	end subroutine

	subroutine iricmi_read_bc_real(type, num, name, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		double precision, intent(out):: val
		integer, intent(out):: ier

		call iricmi_read_bc_real_f2c(type, num, name, val, ier)

	end subroutine

	subroutine iricmi_read_bc_realsingle(type, num, name, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		real, intent(out):: val
		integer, intent(out):: ier

		call iricmi_read_bc_realsingle_f2c(type, num, name, val, ier)

	end subroutine

	subroutine iricmi_read_bc_stringlen(type, num, name, len, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		integer, intent(out):: len
		integer, intent(out):: ier

		call iricmi_read_bc_stringlen_f2c(type, num, name, len, ier)

	end subroutine

	subroutine iricmi_read_bc_string(type, num, name, str, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		character(*), intent(out):: str
		integer, intent(out):: ier

		call iricmi_read_bc_string_f2c(type, num, name, str, ier)

	end subroutine

	subroutine iricmi_read_bc_functional_size(type, num, name, size, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		integer, intent(out):: size
		integer, intent(out):: ier

		call iricmi_read_bc_functional_size_f2c(type, num, name, size, ier)

	end subroutine

	subroutine iricmi_read_bc_functional_vals(type, num, name, x, y, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		double precision, dimension(:), intent(out):: x
		double precision, dimension(:), intent(out):: y
		integer, intent(out):: ier

		call iricmi_read_bc_functional_vals_f2c(type, num, name, x, y, ier)

	end subroutine

	subroutine iricmi_read_bc_functional_valwithname(type, num, name, paramname, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		character(*), intent(in):: paramname
		double precision, dimension(:), intent(out):: val
		integer, intent(out):: ier

		call iricmi_read_bc_functional_valwithname_f2c(type, num, name, paramname, val, ier)

	end subroutine

	subroutine iricmi_read_bc_functional_vals_realsingle(type, num, name, x, y, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		real, dimension(:), intent(out):: x
		real, dimension(:), intent(out):: y
		integer, intent(out):: ier

		call iricmi_read_bc_functional_vals_realsingle_f2c(type, num, name, x, y, ier)

	end subroutine

	subroutine iricmi_read_bc_functional_valwithname_realsingle(type, num, name, paramname, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		character(*), intent(in):: paramname
		real, dimension(:), intent(out):: val
		integer, intent(out):: ier

		call iricmi_read_bc_functional_valwithname_realsingle_f2c(type, num, name, paramname, val, ier)

	end subroutine

	subroutine iricmi_read_complex_count(type, num, ier)
		character(*), intent(in):: type
		integer, intent(out):: num
		integer, intent(out):: ier

		call iricmi_read_complex_count_f2c(type, num, ier)

	end subroutine

	subroutine iricmi_read_complex_integer(type, num, name, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		integer, intent(out):: val
		integer, intent(out):: ier

		call iricmi_read_complex_integer_f2c(type, num, name, val, ier)

	end subroutine

	subroutine iricmi_read_complex_real(type, num, name, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		double precision, intent(out):: val
		integer, intent(out):: ier

		call iricmi_read_complex_real_f2c(type, num, name, val, ier)

	end subroutine

	subroutine iricmi_read_complex_realsingle(type, num, name, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		real, intent(out):: val
		integer, intent(out):: ier

		call iricmi_read_complex_realsingle_f2c(type, num, name, val, ier)

	end subroutine

	subroutine iricmi_read_complex_stringlen(type, num, name, len, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		integer, intent(out):: len
		integer, intent(out):: ier

		call iricmi_read_complex_stringlen_f2c(type, num, name, len, ier)

	end subroutine

	subroutine iricmi_read_complex_string(type, num, name, str, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		character(*), intent(out):: str
		integer, intent(out):: ier

		call iricmi_read_complex_string_f2c(type, num, name, str, ier)

	end subroutine

	subroutine iricmi_read_complex_functional_size(type, num, name, size, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		integer, intent(out):: size
		integer, intent(out):: ier

		call iricmi_read_complex_functional_size_f2c(type, num, name, size, ier)

	end subroutine

	subroutine iricmi_read_complex_functional_vals(type, num, name, x, y, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		double precision, dimension(:), allocatable, intent(out):: x
		double precision, dimension(:), allocatable, intent(out):: y
		integer, intent(out):: ier

		call iricmi_read_complex_functional_vals_f2c(type, num, name, x, y, ier)

	end subroutine

	subroutine iricmi_read_complex_functional_valwithname(type, num, name, paramname, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		character(*), intent(in):: paramname
		double precision, dimension(:), allocatable, intent(out):: val
		integer, intent(out):: ier

		call iricmi_read_complex_functional_valwithname_f2c(type, num, name, paramname, val, ier)

	end subroutine

	subroutine iricmi_read_complex_functional_vals_realsingle(type, num, name, x, y, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		real, dimension(:), allocatable, intent(out):: x
		real, dimension(:), allocatable, intent(out):: y
		integer, intent(out):: ier

		call iricmi_read_complex_functional_vals_realsingle_f2c(type, num, name, x, y, ier)

	end subroutine

	subroutine iricmi_read_complex_functional_valwithname_realsingle(type, num, name, paramname, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: num
		character(*), intent(in):: name
		character(*), intent(in):: paramname
		real, dimension(:), allocatable, intent(out):: val
		integer, intent(out):: ier

		call iricmi_read_complex_functional_valwithname_realsingle_f2c(type, num, name, paramname, val, ier)

	end subroutine

	subroutine iricmi_rin_time(t, ier)
		double precision, intent(out):: t
		integer, intent(out):: ier

		call iricmi_rin_time_f2c(t, ier)

	end subroutine

	subroutine iricmi_rin_dump_interval(interval, ier)
		double precision, intent(out):: interval
		integer, intent(out):: ier

		call iricmi_rin_dump_interval_f2c(interval, ier)

	end subroutine

	subroutine iricmi_rin_integer(name, val, ier)
		character(*), intent(in):: name
		integer, intent(in):: val
		integer, intent(out):: ier

		call iricmi_rin_integer_f2c(name, val, ier)

	end subroutine

	subroutine iricmi_rin_real(name, val, ier)
		character(*), intent(in):: name
		double precision, intent(in):: val
		integer, intent(out):: ier

		call iricmi_rin_real_f2c(name, val, ier)

	end subroutine

	subroutine iricmi_rin_bc_integer(type, index, name, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: index
		character(*), intent(in):: name
		integer, intent(in):: val
		integer, intent(out):: ier

		call iricmi_rin_bc_integer_f2c(type, index, name, val, ier)

	end subroutine

	subroutine iricmi_rin_bc_real(type, index, name, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: index
		character(*), intent(in):: name
		double precision, intent(in):: val
		integer, intent(out):: ier

		call iricmi_rin_bc_real_f2c(type, index, name, val, ier)

	end subroutine

	subroutine iricmi_rin_complex_integer(type, index, name, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: index
		character(*), intent(in):: name
		integer, intent(in):: val
		integer, intent(out):: ier

		call iricmi_rin_complex_integer_f2c(type, index, name, val, ier)

	end subroutine

	subroutine iricmi_rin_complex_real(type, index, name, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: index
		character(*), intent(in):: name
		double precision, intent(in):: val
		integer, intent(out):: ier

		call iricmi_rin_complex_real_f2c(type, index, name, val, ier)

	end subroutine

	subroutine iricmi_rin_grid2d_integer_node(name, vals, ier)
		character(*), intent(in):: name
		integer, dimension(:), allocatable, intent(out):: vals
		integer, intent(out):: ier

		call iricmi_rin_grid2d_integer_node_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_rin_grid2d_real_node(name, vals, ier)
		character(*), intent(in):: name
		double precision, dimension(:), allocatable, intent(out):: vals
		integer, intent(out):: ier

		call iricmi_rin_grid2d_real_node_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_rin_grid2d_integer_cell(name, vals, ier)
		character(*), intent(in):: name
		integer, dimension(:), allocatable, intent(out):: vals
		integer, intent(out):: ier

		call iricmi_rin_grid2d_integer_cell_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_rin_grid2d_real_cell(name, vals, ier)
		character(*), intent(in):: name
		double precision, dimension(:), allocatable, intent(out):: vals
		integer, intent(out):: ier

		call iricmi_rin_grid2d_real_cell_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_rout_time(t, ier)
		double precision, intent(in):: t
		integer, intent(out):: ier

		call iricmi_rout_time_f2c(t, ier)

	end subroutine

	subroutine iricmi_rout_exchange_interval(interval, ier)
		double precision, intent(in):: interval
		integer, intent(out):: ier

		call iricmi_rout_exchange_interval_f2c(interval, ier)

	end subroutine

	subroutine iricmi_rout_dump_interval(interval, ier)
		double precision, intent(in):: interval
		integer, intent(out):: ier

		call iricmi_rout_dump_interval_f2c(interval, ier)

	end subroutine

	subroutine iricmi_rout_integer(name, val, ier)
		character(*), intent(in):: name
		integer, intent(in):: val
		integer, intent(out):: ier

		call iricmi_rout_integer_f2c(name, val, ier)

	end subroutine

	subroutine iricmi_rout_real(name, val, ier)
		character(*), intent(in):: name
		double precision, intent(in):: val
		integer, intent(out):: ier

		call iricmi_rout_real_f2c(name, val, ier)

	end subroutine

	subroutine iricmi_rout_bc_integer(type, index, name, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: index
		character(*), intent(in):: name
		integer, intent(in):: val
		integer, intent(out):: ier

		call iricmi_rout_bc_integer_f2c(type, index, name, val, ier)

	end subroutine

	subroutine iricmi_rout_bc_real(type, index, name, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: index
		character(*), intent(in):: name
		double precision, intent(in):: val
		integer, intent(out):: ier

		call iricmi_rout_bc_real_f2c(type, index, name, val, ier)

	end subroutine

	subroutine iricmi_rout_complex_integer(type, index, name, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: index
		character(*), intent(in):: name
		integer, intent(in):: val
		integer, intent(out):: ier

		call iricmi_rout_complex_integer_f2c(type, index, name, val, ier)

	end subroutine

	subroutine iricmi_rout_complex_real(type, index, name, val, ier)
		character(*), intent(in):: type
		integer, intent(in):: index
		character(*), intent(in):: name
		double precision, intent(in):: val
		integer, intent(out):: ier

		call iricmi_rout_complex_real_f2c(type, index, name, val, ier)

	end subroutine

	subroutine iricmi_rout_grid2d_integer_node_1d(name, vals, ier)
		character(*), intent(in):: name
		integer, dimension(:), intent(in):: vals
		integer, intent(out):: ier

		call iricmi_rout_grid2d_integer_node_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_rout_grid2d_integer_node_2d(name, vals, ier)
		character(*), intent(in):: name
		integer, dimension(:,:), intent(in):: vals
		integer, intent(out):: ier

		call iricmi_rout_grid2d_integer_node_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_rout_grid2d_real_node_1d(name, vals, ier)
		character(*), intent(in):: name
		double precision, dimension(:), intent(in):: vals
		integer, intent(out):: ier

		call iricmi_rout_grid2d_real_node_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_rout_grid2d_real_node_2d(name, vals, ier)
		character(*), intent(in):: name
		double precision, dimension(:,:), intent(in):: vals
		integer, intent(out):: ier

		call iricmi_rout_grid2d_real_node_f2c(name, vals, ier)

	end subroutine


	subroutine iricmi_rout_grid2d_integer_cell_1d(name, vals, ier)
		character(*), intent(in):: name
		integer, dimension(:), intent(in):: vals
		integer, intent(out):: ier

		call iricmi_rout_grid2d_integer_cell_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_rout_grid2d_integer_cell_2d(name, vals, ier)
		character(*), intent(in):: name
		integer, dimension(:,:), intent(in):: vals
		integer, intent(out):: ier

		call iricmi_rout_grid2d_integer_cell_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_rout_grid2d_real_cell_1d(name, vals, ier)
		character(*), intent(in):: name
		double precision, dimension(:), intent(in):: vals
		integer, intent(out):: ier

		call iricmi_rout_grid2d_real_cell_f2c(name, vals, ier)

	end subroutine

	subroutine iricmi_rout_grid2d_real_cell_2d(name, vals, ier)
		character(*), intent(in):: name
		double precision, dimension(:,:), intent(in):: vals
		integer, intent(out):: ier

		call iricmi_rout_grid2d_real_cell_f2c(name, vals, ier)

	end subroutine

end module
