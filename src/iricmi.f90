module iricmi
  implicit none

  integer, parameter:: IRICMI_GEO_UNKNOWN = 0
  integer, parameter:: IRICMI_GEO_POLYGON = 1
  integer, parameter:: IRICMI_GEO_RIVERSURVEY = 2
  integer, parameter:: IRICMI_GEO_POINTMAP = 3
  integer, parameter:: IRICMI_GEO_POLYLINE = 4
  interface iricmi_read_bc_indices_withgridid
    module procedure iricmi_read_bc_indices_withgridid_1d
    module procedure iricmi_read_bc_indices_withgridid_2d
    module procedure iricmi_read_bc_indices_withgridid_3d
  end interface

  interface iricmi_read_bc_functional_withgridid
    module procedure iricmi_read_bc_functional_withgridid_1d
    module procedure iricmi_read_bc_functional_withgridid_2d
    module procedure iricmi_read_bc_functional_withgridid_3d
  end interface

  interface iricmi_read_bc_functionalwithname_withgridid
    module procedure iricmi_read_bc_functionalwithname_withgridid_1d
    module procedure iricmi_read_bc_functionalwithname_withgridid_2d
    module procedure iricmi_read_bc_functionalwithname_withgridid_3d
  end interface

  interface iricmi_read_bc_functional_realsingle_withgridid
    module procedure iricmi_read_bc_functional_realsingle_withgridid_1d
    module procedure iricmi_read_bc_functional_realsingle_withgridid_2d
    module procedure iricmi_read_bc_functional_realsingle_withgridid_3d
  end interface

  interface iricmi_read_bc_functionalwithname_realsingle_withgridid
    module procedure iricmi_read_bc_functionalwithname_realsingle_withgridid_1d
    module procedure iricmi_read_bc_functionalwithname_realsingle_withgridid_2d
    module procedure iricmi_read_bc_functionalwithname_realsingle_withgridid_3d
  end interface

  interface iricmi_read_functional
    module procedure iricmi_read_functional_1d
    module procedure iricmi_read_functional_2d
    module procedure iricmi_read_functional_3d
  end interface

  interface iricmi_read_functionalwithname
    module procedure iricmi_read_functionalwithname_1d
    module procedure iricmi_read_functionalwithname_2d
    module procedure iricmi_read_functionalwithname_3d
  end interface

  interface iricmi_read_functional_realsingle
    module procedure iricmi_read_functional_realsingle_1d
    module procedure iricmi_read_functional_realsingle_2d
    module procedure iricmi_read_functional_realsingle_3d
  end interface

  interface iricmi_read_functionalwithname_realsingle
    module procedure iricmi_read_functionalwithname_realsingle_1d
    module procedure iricmi_read_functionalwithname_realsingle_2d
    module procedure iricmi_read_functionalwithname_realsingle_3d
  end interface

  interface iricmi_read_grid_complex_node_withgridid
    module procedure iricmi_read_grid_complex_node_withgridid_1d
    module procedure iricmi_read_grid_complex_node_withgridid_2d
    module procedure iricmi_read_grid_complex_node_withgridid_3d
  end interface

  interface iricmi_read_grid_complex_cell_withgridid
    module procedure iricmi_read_grid_complex_cell_withgridid_1d
    module procedure iricmi_read_grid_complex_cell_withgridid_2d
    module procedure iricmi_read_grid_complex_cell_withgridid_3d
  end interface

  interface iricmi_read_grid2d_coords_withgridid
    module procedure iricmi_read_grid2d_coords_withgridid_1d
    module procedure iricmi_read_grid2d_coords_withgridid_2d
  end interface

  interface iricmi_read_grid3d_coords_withgridid
    module procedure iricmi_read_grid3d_coords_withgridid_1d
    module procedure iricmi_read_grid3d_coords_withgridid_3d
  end interface

  interface iricmi_read_grid_real_node_withgridid
    module procedure iricmi_read_grid_real_node_withgridid_1d
    module procedure iricmi_read_grid_real_node_withgridid_2d
    module procedure iricmi_read_grid_real_node_withgridid_3d
  end interface

  interface iricmi_read_grid_integer_node_withgridid
    module procedure iricmi_read_grid_integer_node_withgridid_1d
    module procedure iricmi_read_grid_integer_node_withgridid_2d
    module procedure iricmi_read_grid_integer_node_withgridid_3d
  end interface

  interface iricmi_read_grid_real_cell_withgridid
    module procedure iricmi_read_grid_real_cell_withgridid_1d
    module procedure iricmi_read_grid_real_cell_withgridid_2d
    module procedure iricmi_read_grid_real_cell_withgridid_3d
  end interface

  interface iricmi_read_grid_integer_cell_withgridid
    module procedure iricmi_read_grid_integer_cell_withgridid_1d
    module procedure iricmi_read_grid_integer_cell_withgridid_2d
    module procedure iricmi_read_grid_integer_cell_withgridid_3d
  end interface

  interface iricmi_read_grid_functional_integer_node_withgridid
    module procedure iricmi_read_grid_functional_integer_node_withgridid_1d
    module procedure iricmi_read_grid_functional_integer_node_withgridid_2d
    module procedure iricmi_read_grid_functional_integer_node_withgridid_3d
  end interface

  interface iricmi_read_grid_functional_real_node_withgridid
    module procedure iricmi_read_grid_functional_real_node_withgridid_1d
    module procedure iricmi_read_grid_functional_real_node_withgridid_2d
    module procedure iricmi_read_grid_functional_real_node_withgridid_3d
  end interface

  interface iricmi_read_grid_functional_integer_cell_withgridid
    module procedure iricmi_read_grid_functional_integer_cell_withgridid_1d
    module procedure iricmi_read_grid_functional_integer_cell_withgridid_2d
    module procedure iricmi_read_grid_functional_integer_cell_withgridid_3d
  end interface

  interface iricmi_read_grid_functional_real_cell_withgridid
    module procedure iricmi_read_grid_functional_real_cell_withgridid_1d
    module procedure iricmi_read_grid_functional_real_cell_withgridid_2d
    module procedure iricmi_read_grid_functional_real_cell_withgridid_3d
  end interface

  interface iricmi_write_grid2d_coords_withgridid
    module procedure iricmi_write_grid2d_coords_withgridid_1d
    module procedure iricmi_write_grid2d_coords_withgridid_2d
  end interface

  interface iricmi_write_grid3d_coords_withgridid
    module procedure iricmi_write_grid3d_coords_withgridid_1d
    module procedure iricmi_write_grid3d_coords_withgridid_3d
  end interface

  interface iricmi_write_namedgrid2d_coords_withgridid
    module procedure iricmi_write_namedgrid2d_coords_withgridid_1d
    module procedure iricmi_write_namedgrid2d_coords_withgridid_2d
  end interface

  interface iricmi_write_namedgrid3d_coords_withgridid
    module procedure iricmi_write_namedgrid3d_coords_withgridid_1d
    module procedure iricmi_write_namedgrid3d_coords_withgridid_3d
  end interface

  interface iricmi_read_complex_functional
    module procedure iricmi_read_complex_functional_1d
    module procedure iricmi_read_complex_functional_2d
    module procedure iricmi_read_complex_functional_3d
  end interface

  interface iricmi_read_complex_functionalwithname
    module procedure iricmi_read_complex_functionalwithname_1d
    module procedure iricmi_read_complex_functionalwithname_2d
    module procedure iricmi_read_complex_functionalwithname_3d
  end interface

  interface iricmi_read_complex_functional_realsingle
    module procedure iricmi_read_complex_functional_realsingle_1d
    module procedure iricmi_read_complex_functional_realsingle_2d
    module procedure iricmi_read_complex_functional_realsingle_3d
  end interface

  interface iricmi_read_complex_functionalwithname_realsingle
    module procedure iricmi_read_complex_functionalwithname_realsingle_1d
    module procedure iricmi_read_complex_functionalwithname_realsingle_2d
    module procedure iricmi_read_complex_functionalwithname_realsingle_3d
  end interface

  interface iricmi_read_bc_indices
    module procedure iricmi_read_bc_indices_1d
    module procedure iricmi_read_bc_indices_2d
    module procedure iricmi_read_bc_indices_3d
  end interface

  interface iricmi_read_bc_functional
    module procedure iricmi_read_bc_functional_1d
    module procedure iricmi_read_bc_functional_2d
    module procedure iricmi_read_bc_functional_3d
  end interface

  interface iricmi_read_bc_functionalwithname
    module procedure iricmi_read_bc_functionalwithname_1d
    module procedure iricmi_read_bc_functionalwithname_2d
    module procedure iricmi_read_bc_functionalwithname_3d
  end interface

  interface iricmi_read_bc_functional_realsingle
    module procedure iricmi_read_bc_functional_realsingle_1d
    module procedure iricmi_read_bc_functional_realsingle_2d
    module procedure iricmi_read_bc_functional_realsingle_3d
  end interface

  interface iricmi_read_bc_functionalwithname_realsingle
    module procedure iricmi_read_bc_functionalwithname_realsingle_1d
    module procedure iricmi_read_bc_functionalwithname_realsingle_2d
    module procedure iricmi_read_bc_functionalwithname_realsingle_3d
  end interface

  interface iricmi_read_grid_complex_node
    module procedure iricmi_read_grid_complex_node_1d
    module procedure iricmi_read_grid_complex_node_2d
    module procedure iricmi_read_grid_complex_node_3d
  end interface

  interface iricmi_read_grid_complex_cell
    module procedure iricmi_read_grid_complex_cell_1d
    module procedure iricmi_read_grid_complex_cell_2d
    module procedure iricmi_read_grid_complex_cell_3d
  end interface

  interface iricmi_read_grid2d_coords
    module procedure iricmi_read_grid2d_coords_1d
    module procedure iricmi_read_grid2d_coords_2d
  end interface

  interface iricmi_read_grid3d_coords
    module procedure iricmi_read_grid3d_coords_1d
    module procedure iricmi_read_grid3d_coords_3d
  end interface

  interface iricmi_read_grid_real_node
    module procedure iricmi_read_grid_real_node_1d
    module procedure iricmi_read_grid_real_node_2d
    module procedure iricmi_read_grid_real_node_3d
  end interface

  interface iricmi_read_grid_integer_node
    module procedure iricmi_read_grid_integer_node_1d
    module procedure iricmi_read_grid_integer_node_2d
    module procedure iricmi_read_grid_integer_node_3d
  end interface

  interface iricmi_read_grid_real_cell
    module procedure iricmi_read_grid_real_cell_1d
    module procedure iricmi_read_grid_real_cell_2d
    module procedure iricmi_read_grid_real_cell_3d
  end interface

  interface iricmi_read_grid_integer_cell
    module procedure iricmi_read_grid_integer_cell_1d
    module procedure iricmi_read_grid_integer_cell_2d
    module procedure iricmi_read_grid_integer_cell_3d
  end interface

  interface iricmi_read_grid_functional_integer_node
    module procedure iricmi_read_grid_functional_integer_node_1d
    module procedure iricmi_read_grid_functional_integer_node_2d
    module procedure iricmi_read_grid_functional_integer_node_3d
  end interface

  interface iricmi_read_grid_functional_real_node
    module procedure iricmi_read_grid_functional_real_node_1d
    module procedure iricmi_read_grid_functional_real_node_2d
    module procedure iricmi_read_grid_functional_real_node_3d
  end interface

  interface iricmi_read_grid_functional_integer_cell
    module procedure iricmi_read_grid_functional_integer_cell_1d
    module procedure iricmi_read_grid_functional_integer_cell_2d
    module procedure iricmi_read_grid_functional_integer_cell_3d
  end interface

  interface iricmi_read_grid_functional_real_cell
    module procedure iricmi_read_grid_functional_real_cell_1d
    module procedure iricmi_read_grid_functional_real_cell_2d
    module procedure iricmi_read_grid_functional_real_cell_3d
  end interface

  interface iricmi_write_grid2d_coords
    module procedure iricmi_write_grid2d_coords_1d
    module procedure iricmi_write_grid2d_coords_2d
  end interface

  interface iricmi_write_grid3d_coords
    module procedure iricmi_write_grid3d_coords_1d
    module procedure iricmi_write_grid3d_coords_3d
  end interface

  interface iricmi_write_namedgrid2d_coords
    module procedure iricmi_write_namedgrid2d_coords_1d
    module procedure iricmi_write_namedgrid2d_coords_2d
  end interface

  interface iricmi_write_namedgrid3d_coords
    module procedure iricmi_write_namedgrid3d_coords_1d
    module procedure iricmi_write_namedgrid3d_coords_3d
  end interface

  interface iricmi_rin_grid_cell_integer
    module procedure iricmi_rin_grid_cell_integer_1d
    module procedure iricmi_rin_grid_cell_integer_2d
    module procedure iricmi_rin_grid_cell_integer_3d
  end interface

  interface iricmi_rin_grid_cell_real
    module procedure iricmi_rin_grid_cell_real_1d
    module procedure iricmi_rin_grid_cell_real_2d
    module procedure iricmi_rin_grid_cell_real_3d
  end interface

  interface iricmi_rout_grid_cell_integer
    module procedure iricmi_rout_grid_cell_integer_1d
    module procedure iricmi_rout_grid_cell_integer_2d
    module procedure iricmi_rout_grid_cell_integer_3d
  end interface

  interface iricmi_rout_grid_cell_real
    module procedure iricmi_rout_grid_cell_real_1d
    module procedure iricmi_rout_grid_cell_real_2d
    module procedure iricmi_rout_grid_cell_real_3d
  end interface

  interface iricmi_rin_grid2d_coords
    module procedure iricmi_rin_grid2d_coords_1d
    module procedure iricmi_rin_grid2d_coords_2d
  end interface

  interface iricmi_rin_grid3d_coords
    module procedure iricmi_rin_grid3d_coords_1d
    module procedure iricmi_rin_grid3d_coords_3d
  end interface

  interface iricmi_rout_grid2d_coords
    module procedure iricmi_rout_grid2d_coords_1d
    module procedure iricmi_rout_grid2d_coords_2d
  end interface

  interface iricmi_rout_grid3d_coords
    module procedure iricmi_rout_grid3d_coords_1d
    module procedure iricmi_rout_grid3d_coords_3d
  end interface

  interface iricmi_rin_grid_iface_integer
    module procedure iricmi_rin_grid_iface_integer_1d
    module procedure iricmi_rin_grid_iface_integer_2d
    module procedure iricmi_rin_grid_iface_integer_3d
  end interface

  interface iricmi_rin_grid_iface_real
    module procedure iricmi_rin_grid_iface_real_1d
    module procedure iricmi_rin_grid_iface_real_2d
    module procedure iricmi_rin_grid_iface_real_3d
  end interface

  interface iricmi_rout_grid_iface_integer
    module procedure iricmi_rout_grid_iface_integer_1d
    module procedure iricmi_rout_grid_iface_integer_2d
    module procedure iricmi_rout_grid_iface_integer_3d
  end interface

  interface iricmi_rout_grid_iface_real
    module procedure iricmi_rout_grid_iface_real_1d
    module procedure iricmi_rout_grid_iface_real_2d
    module procedure iricmi_rout_grid_iface_real_3d
  end interface

  interface iricmi_rin_grid_jface_integer
    module procedure iricmi_rin_grid_jface_integer_1d
    module procedure iricmi_rin_grid_jface_integer_2d
    module procedure iricmi_rin_grid_jface_integer_3d
  end interface

  interface iricmi_rin_grid_jface_real
    module procedure iricmi_rin_grid_jface_real_1d
    module procedure iricmi_rin_grid_jface_real_2d
    module procedure iricmi_rin_grid_jface_real_3d
  end interface

  interface iricmi_rout_grid_jface_integer
    module procedure iricmi_rout_grid_jface_integer_1d
    module procedure iricmi_rout_grid_jface_integer_2d
    module procedure iricmi_rout_grid_jface_integer_3d
  end interface

  interface iricmi_rout_grid_jface_real
    module procedure iricmi_rout_grid_jface_real_1d
    module procedure iricmi_rout_grid_jface_real_2d
    module procedure iricmi_rout_grid_jface_real_3d
  end interface

  interface iricmi_rin_grid_kface_integer
    module procedure iricmi_rin_grid_kface_integer_1d
    module procedure iricmi_rin_grid_kface_integer_2d
    module procedure iricmi_rin_grid_kface_integer_3d
  end interface

  interface iricmi_rin_grid_kface_real
    module procedure iricmi_rin_grid_kface_real_1d
    module procedure iricmi_rin_grid_kface_real_2d
    module procedure iricmi_rin_grid_kface_real_3d
  end interface

  interface iricmi_rout_grid_kface_integer
    module procedure iricmi_rout_grid_kface_integer_1d
    module procedure iricmi_rout_grid_kface_integer_2d
    module procedure iricmi_rout_grid_kface_integer_3d
  end interface

  interface iricmi_rout_grid_kface_real
    module procedure iricmi_rout_grid_kface_real_1d
    module procedure iricmi_rout_grid_kface_real_2d
    module procedure iricmi_rout_grid_kface_real_3d
  end interface

  interface iricmi_rin_grid_node_integer
    module procedure iricmi_rin_grid_node_integer_1d
    module procedure iricmi_rin_grid_node_integer_2d
    module procedure iricmi_rin_grid_node_integer_3d
  end interface

  interface iricmi_rin_grid_node_real
    module procedure iricmi_rin_grid_node_real_1d
    module procedure iricmi_rin_grid_node_real_2d
    module procedure iricmi_rin_grid_node_real_3d
  end interface

  interface iricmi_rout_grid_node_integer
    module procedure iricmi_rout_grid_node_integer_1d
    module procedure iricmi_rout_grid_node_integer_2d
    module procedure iricmi_rout_grid_node_integer_3d
  end interface

  interface iricmi_rout_grid_node_real
    module procedure iricmi_rout_grid_node_real_1d
    module procedure iricmi_rout_grid_node_real_2d
    module procedure iricmi_rout_grid_node_real_3d
  end interface

  interface iricmi_rin_grid_cell_integer_withgridid
    module procedure iricmi_rin_grid_cell_integer_withgridid_1d
    module procedure iricmi_rin_grid_cell_integer_withgridid_2d
    module procedure iricmi_rin_grid_cell_integer_withgridid_3d
  end interface

  interface iricmi_rin_grid_cell_real_withgridid
    module procedure iricmi_rin_grid_cell_real_withgridid_1d
    module procedure iricmi_rin_grid_cell_real_withgridid_2d
    module procedure iricmi_rin_grid_cell_real_withgridid_3d
  end interface

  interface iricmi_rout_grid_cell_integer_withgridid
    module procedure iricmi_rout_grid_cell_integer_withgridid_1d
    module procedure iricmi_rout_grid_cell_integer_withgridid_2d
    module procedure iricmi_rout_grid_cell_integer_withgridid_3d
  end interface

  interface iricmi_rout_grid_cell_real_withgridid
    module procedure iricmi_rout_grid_cell_real_withgridid_1d
    module procedure iricmi_rout_grid_cell_real_withgridid_2d
    module procedure iricmi_rout_grid_cell_real_withgridid_3d
  end interface

  interface iricmi_rin_grid2d_coords_withgridid
    module procedure iricmi_rin_grid2d_coords_withgridid_1d
    module procedure iricmi_rin_grid2d_coords_withgridid_2d
  end interface

  interface iricmi_rin_grid3d_coords_withgridid
    module procedure iricmi_rin_grid3d_coords_withgridid_1d
    module procedure iricmi_rin_grid3d_coords_withgridid_3d
  end interface

  interface iricmi_rout_grid2d_coords_withgridid
    module procedure iricmi_rout_grid2d_coords_withgridid_1d
    module procedure iricmi_rout_grid2d_coords_withgridid_2d
  end interface

  interface iricmi_rout_grid3d_coords_withgridid
    module procedure iricmi_rout_grid3d_coords_withgridid_1d
    module procedure iricmi_rout_grid3d_coords_withgridid_3d
  end interface

  interface iricmi_rin_grid_iface_integer_withgridid
    module procedure iricmi_rin_grid_iface_integer_withgridid_1d
    module procedure iricmi_rin_grid_iface_integer_withgridid_2d
    module procedure iricmi_rin_grid_iface_integer_withgridid_3d
  end interface

  interface iricmi_rin_grid_iface_real_withgridid
    module procedure iricmi_rin_grid_iface_real_withgridid_1d
    module procedure iricmi_rin_grid_iface_real_withgridid_2d
    module procedure iricmi_rin_grid_iface_real_withgridid_3d
  end interface

  interface iricmi_rout_grid_iface_integer_withgridid
    module procedure iricmi_rout_grid_iface_integer_withgridid_1d
    module procedure iricmi_rout_grid_iface_integer_withgridid_2d
    module procedure iricmi_rout_grid_iface_integer_withgridid_3d
  end interface

  interface iricmi_rout_grid_iface_real_withgridid
    module procedure iricmi_rout_grid_iface_real_withgridid_1d
    module procedure iricmi_rout_grid_iface_real_withgridid_2d
    module procedure iricmi_rout_grid_iface_real_withgridid_3d
  end interface

  interface iricmi_rin_grid_jface_integer_withgridid
    module procedure iricmi_rin_grid_jface_integer_withgridid_1d
    module procedure iricmi_rin_grid_jface_integer_withgridid_2d
    module procedure iricmi_rin_grid_jface_integer_withgridid_3d
  end interface

  interface iricmi_rin_grid_jface_real_withgridid
    module procedure iricmi_rin_grid_jface_real_withgridid_1d
    module procedure iricmi_rin_grid_jface_real_withgridid_2d
    module procedure iricmi_rin_grid_jface_real_withgridid_3d
  end interface

  interface iricmi_rout_grid_jface_integer_withgridid
    module procedure iricmi_rout_grid_jface_integer_withgridid_1d
    module procedure iricmi_rout_grid_jface_integer_withgridid_2d
    module procedure iricmi_rout_grid_jface_integer_withgridid_3d
  end interface

  interface iricmi_rout_grid_jface_real_withgridid
    module procedure iricmi_rout_grid_jface_real_withgridid_1d
    module procedure iricmi_rout_grid_jface_real_withgridid_2d
    module procedure iricmi_rout_grid_jface_real_withgridid_3d
  end interface

  interface iricmi_rin_grid_kface_integer_withgridid
    module procedure iricmi_rin_grid_kface_integer_withgridid_1d
    module procedure iricmi_rin_grid_kface_integer_withgridid_2d
    module procedure iricmi_rin_grid_kface_integer_withgridid_3d
  end interface

  interface iricmi_rin_grid_kface_real_withgridid
    module procedure iricmi_rin_grid_kface_real_withgridid_1d
    module procedure iricmi_rin_grid_kface_real_withgridid_2d
    module procedure iricmi_rin_grid_kface_real_withgridid_3d
  end interface

  interface iricmi_rout_grid_kface_integer_withgridid
    module procedure iricmi_rout_grid_kface_integer_withgridid_1d
    module procedure iricmi_rout_grid_kface_integer_withgridid_2d
    module procedure iricmi_rout_grid_kface_integer_withgridid_3d
  end interface

  interface iricmi_rout_grid_kface_real_withgridid
    module procedure iricmi_rout_grid_kface_real_withgridid_1d
    module procedure iricmi_rout_grid_kface_real_withgridid_2d
    module procedure iricmi_rout_grid_kface_real_withgridid_3d
  end interface

  interface iricmi_rin_grid_node_integer_withgridid
    module procedure iricmi_rin_grid_node_integer_withgridid_1d
    module procedure iricmi_rin_grid_node_integer_withgridid_2d
    module procedure iricmi_rin_grid_node_integer_withgridid_3d
  end interface

  interface iricmi_rin_grid_node_real_withgridid
    module procedure iricmi_rin_grid_node_real_withgridid_1d
    module procedure iricmi_rin_grid_node_real_withgridid_2d
    module procedure iricmi_rin_grid_node_real_withgridid_3d
  end interface

  interface iricmi_rout_grid_node_integer_withgridid
    module procedure iricmi_rout_grid_node_integer_withgridid_1d
    module procedure iricmi_rout_grid_node_integer_withgridid_2d
    module procedure iricmi_rout_grid_node_integer_withgridid_3d
  end interface

  interface iricmi_rout_grid_node_real_withgridid
    module procedure iricmi_rout_grid_node_real_withgridid_1d
    module procedure iricmi_rout_grid_node_real_withgridid_2d
    module procedure iricmi_rout_grid_node_real_withgridid_3d
  end interface


contains


  ! from iricmi_bc.h

  subroutine iricmi_read_bc_count_withgridid(gid, type, num, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(out):: num
    integer, intent(out):: ier

    call iricmi_read_bc_count_withgridid_f2c &
      (gid, type, num, ier)

  end subroutine

  subroutine iricmi_read_bc_indicessize_withgridid(gid, type, num, size, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    integer, intent(out):: size
    integer, intent(out):: ier

    call iricmi_read_bc_indicessize_withgridid_f2c &
      (gid, type, num, size, ier)

  end subroutine

  subroutine iricmi_read_bc_indices_withgridid_1d(gid, type, num, idx_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    integer, dimension(:), intent(out):: idx_arr
    integer, intent(out):: ier

    call iricmi_read_bc_indices_withgridid_f2c &
      (gid, type, num, idx_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_indices_withgridid_2d(gid, type, num, idx_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    integer, dimension(:,:), intent(out):: idx_arr
    integer, intent(out):: ier

    call iricmi_read_bc_indices_withgridid_f2c &
      (gid, type, num, idx_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_indices_withgridid_3d(gid, type, num, idx_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    integer, dimension(:,:,:), intent(out):: idx_arr
    integer, intent(out):: ier

    call iricmi_read_bc_indices_withgridid_f2c &
      (gid, type, num, idx_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_integer_withgridid(gid, type, num, name, value, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    integer, intent(out):: value
    integer, intent(out):: ier

    call iricmi_read_bc_integer_withgridid_f2c &
      (gid, type, num, name, value, ier)

  end subroutine

  subroutine iricmi_read_bc_real_withgridid(gid, type, num, name, value, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, intent(out):: value
    integer, intent(out):: ier

    call iricmi_read_bc_real_withgridid_f2c &
      (gid, type, num, name, value, ier)

  end subroutine

  subroutine iricmi_read_bc_realsingle_withgridid(gid, type, num, name, value, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    real, intent(out):: value
    integer, intent(out):: ier

    call iricmi_read_bc_realsingle_withgridid_f2c &
      (gid, type, num, name, value, ier)

  end subroutine

  subroutine iricmi_read_bc_stringlen_withgridid(gid, type, num, name, length, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    integer, intent(out):: length
    integer, intent(out):: ier

    call iricmi_read_bc_stringlen_withgridid_f2c &
      (gid, type, num, name, length, ier)

  end subroutine

  subroutine iricmi_read_bc_string_withgridid(gid, type, num, name, strvalue, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(out):: strvalue
    integer, intent(out):: ier

    call iricmi_read_bc_string_withgridid_f2c &
      (gid, type, num, name, strvalue, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalsize_withgridid(gid, type, num, name, size, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    integer, intent(out):: size
    integer, intent(out):: ier

    call iricmi_read_bc_functionalsize_withgridid_f2c &
      (gid, type, num, name, size, ier)

  end subroutine

  subroutine iricmi_read_bc_functional_withgridid_1d(gid, type, num, name, x_arr, y_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, dimension(:), intent(out):: x_arr
    double precision, dimension(:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functional_withgridid_f2c &
      (gid, type, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functional_withgridid_2d(gid, type, num, name, x_arr, y_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(out):: x_arr
    double precision, dimension(:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functional_withgridid_f2c &
      (gid, type, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functional_withgridid_3d(gid, type, num, name, x_arr, y_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(out):: x_arr
    double precision, dimension(:,:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functional_withgridid_f2c &
      (gid, type, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functional_withgridid_4d(gid, type, num, name, x_arr, y_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, dimension(:,:,:,:), intent(out):: x_arr
    double precision, dimension(:,:,:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functional_withgridid_f2c &
      (gid, type, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_withgridid_1d(gid, type, num, name, paramname, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    double precision, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_withgridid_f2c &
      (gid, type, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_withgridid_2d(gid, type, num, name, paramname, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    double precision, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_withgridid_f2c &
      (gid, type, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_withgridid_3d(gid, type, num, name, paramname, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    double precision, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_withgridid_f2c &
      (gid, type, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_withgridid_4d(gid, type, num, name, paramname, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    double precision, dimension(:,:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_withgridid_f2c &
      (gid, type, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functional_realsingle_withgridid_1d(gid, type, num, name, x_arr, y_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    real, dimension(:), intent(out):: x_arr
    real, dimension(:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functional_realsingle_withgridid_f2c &
      (gid, type, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functional_realsingle_withgridid_2d(gid, type, num, name, x_arr, y_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    real, dimension(:,:), intent(out):: x_arr
    real, dimension(:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functional_realsingle_withgridid_f2c &
      (gid, type, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functional_realsingle_withgridid_3d(gid, type, num, name, x_arr, y_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    real, dimension(:,:,:), intent(out):: x_arr
    real, dimension(:,:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functional_realsingle_withgridid_f2c &
      (gid, type, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functional_realsingle_withgridid_4d(gid, type, num, name, x_arr, y_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    real, dimension(:,:,:,:), intent(out):: x_arr
    real, dimension(:,:,:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functional_realsingle_withgridid_f2c &
      (gid, type, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_realsingle_withgridid_1d(gid, type, num, name, paramname, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    real, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_realsingle_withgridid_f2c &
      (gid, type, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_realsingle_withgridid_2d(gid, type, num, name, paramname, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    real, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_realsingle_withgridid_f2c &
      (gid, type, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_realsingle_withgridid_3d(gid, type, num, name, paramname, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    real, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_realsingle_withgridid_f2c &
      (gid, type, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_realsingle_withgridid_4d(gid, type, num, name, paramname, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    real, dimension(:,:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_realsingle_withgridid_f2c &
      (gid, type, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_stringlen_withgridid(gid, type, num, name, paramname, length, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    integer, intent(out):: length
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_stringlen_withgridid_f2c &
      (gid, type, num, name, paramname, length, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_string_withgridid(gid, type, num, name, paramname, strvalue, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    character(*), intent(out):: strvalue
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_string_withgridid_f2c &
      (gid, type, num, name, paramname, strvalue, ier)

  end subroutine

  subroutine iricmi_rin_bc_integer_withgridid(gid, type, num, name, val, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    integer, intent(in):: val
    integer, intent(out):: ier

    call iricmi_rin_bc_integer_withgridid_f2c &
      (gid, type, num, name, val, ier)

  end subroutine

  subroutine iricmi_rin_bc_real_withgridid(gid, type, num, name, val, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, intent(in):: val
    integer, intent(out):: ier

    call iricmi_rin_bc_real_withgridid_f2c &
      (gid, type, num, name, val, ier)

  end subroutine

  subroutine iricmi_rout_bc_integer_withgridid(gid, type, num, name, val, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    integer, intent(in):: val
    integer, intent(out):: ier

    call iricmi_rout_bc_integer_withgridid_f2c &
      (gid, type, num, name, val, ier)

  end subroutine

  subroutine iricmi_rout_bc_real_withgridid(gid, type, num, name, val, ier)
    integer, intent(in):: gid
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, intent(in):: val
    integer, intent(out):: ier

    call iricmi_rout_bc_real_withgridid_f2c &
      (gid, type, num, name, val, ier)

  end subroutine



  ! from iricmi_calc.h

  subroutine iricmi_calc_init(ier)
    integer, intent(out):: ier

    call iricmi_calc_init_f2c &
      (ier)

  end subroutine

  subroutine iricmi_calc_terminate(ier)
    integer, intent(out):: ier

    call iricmi_calc_terminate_f2c &
      (ier)

  end subroutine

  subroutine iricmi_calc_sync_receive(ier)
    integer, intent(out):: ier

    call iricmi_calc_sync_receive_f2c &
      (ier)

  end subroutine

  subroutine iricmi_calc_sync_send(ier)
    integer, intent(out):: ier

    call iricmi_calc_sync_send_f2c &
      (ier)

  end subroutine

  subroutine iricmi_calc_dump(ier)
    integer, intent(out):: ier

    call iricmi_calc_dump_f2c &
      (ier)

  end subroutine



  ! from iricmi_cc.h

  subroutine iricmi_read_integer(name, value, ier)
    character(*), intent(in):: name
    integer, intent(out):: value
    integer, intent(out):: ier

    call iricmi_read_integer_f2c &
      (name, value, ier)

  end subroutine

  subroutine iricmi_read_real(name, value, ier)
    character(*), intent(in):: name
    double precision, intent(out):: value
    integer, intent(out):: ier

    call iricmi_read_real_f2c &
      (name, value, ier)

  end subroutine

  subroutine iricmi_read_realsingle(name, value, ier)
    character(*), intent(in):: name
    real, intent(out):: value
    integer, intent(out):: ier

    call iricmi_read_realsingle_f2c &
      (name, value, ier)

  end subroutine

  subroutine iricmi_read_stringlen(name, length, ier)
    character(*), intent(in):: name
    integer, intent(out):: length
    integer, intent(out):: ier

    call iricmi_read_stringlen_f2c &
      (name, length, ier)

  end subroutine

  subroutine iricmi_read_string(name, strvalue, ier)
    character(*), intent(in):: name
    character(*), intent(out):: strvalue
    integer, intent(out):: ier

    call iricmi_read_string_f2c &
      (name, strvalue, ier)

  end subroutine

  subroutine iricmi_read_functionalsize(name, size, ier)
    character(*), intent(in):: name
    integer, intent(out):: size
    integer, intent(out):: ier

    call iricmi_read_functionalsize_f2c &
      (name, size, ier)

  end subroutine

  subroutine iricmi_read_functional_1d(name, x_arr, y_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:), intent(out):: x_arr
    double precision, dimension(:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_functional_f2c &
      (name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_functional_2d(name, x_arr, y_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(out):: x_arr
    double precision, dimension(:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_functional_f2c &
      (name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_functional_3d(name, x_arr, y_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(out):: x_arr
    double precision, dimension(:,:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_functional_f2c &
      (name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_functional_4d(name, x_arr, y_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:,:,:), intent(out):: x_arr
    double precision, dimension(:,:,:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_functional_f2c &
      (name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_functionalwithname_1d(name, paramname, v_arr, ier)
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    double precision, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_functionalwithname_f2c &
      (name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_functionalwithname_2d(name, paramname, v_arr, ier)
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    double precision, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_functionalwithname_f2c &
      (name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_functionalwithname_3d(name, paramname, v_arr, ier)
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    double precision, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_functionalwithname_f2c &
      (name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_functionalwithname_4d(name, paramname, v_arr, ier)
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    double precision, dimension(:,:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_functionalwithname_f2c &
      (name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_functional_realsingle_1d(name, x_arr, y_arr, ier)
    character(*), intent(in):: name
    real, dimension(:), intent(out):: x_arr
    real, dimension(:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_functional_realsingle_f2c &
      (name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_functional_realsingle_2d(name, x_arr, y_arr, ier)
    character(*), intent(in):: name
    real, dimension(:,:), intent(out):: x_arr
    real, dimension(:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_functional_realsingle_f2c &
      (name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_functional_realsingle_3d(name, x_arr, y_arr, ier)
    character(*), intent(in):: name
    real, dimension(:,:,:), intent(out):: x_arr
    real, dimension(:,:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_functional_realsingle_f2c &
      (name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_functional_realsingle_4d(name, x_arr, y_arr, ier)
    character(*), intent(in):: name
    real, dimension(:,:,:,:), intent(out):: x_arr
    real, dimension(:,:,:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_functional_realsingle_f2c &
      (name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_functionalwithname_realsingle_1d(name, paramname, v_arr, ier)
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    real, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_functionalwithname_realsingle_f2c &
      (name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_functionalwithname_realsingle_2d(name, paramname, v_arr, ier)
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    real, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_functionalwithname_realsingle_f2c &
      (name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_functionalwithname_realsingle_3d(name, paramname, v_arr, ier)
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    real, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_functionalwithname_realsingle_f2c &
      (name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_functionalwithname_realsingle_4d(name, paramname, v_arr, ier)
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    real, dimension(:,:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_functionalwithname_realsingle_f2c &
      (name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_functionalwithname_string(name, paramname, strvalue, ier)
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    character(*), intent(out):: strvalue
    integer, intent(out):: ier

    call iricmi_read_functionalwithname_string_f2c &
      (name, paramname, strvalue, ier)

  end subroutine

  subroutine iricmi_read_functionalwithname_stringlen(name, paramname, length, ier)
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    integer, intent(out):: length
    integer, intent(out):: ier

    call iricmi_read_functionalwithname_stringlen_f2c &
      (name, paramname, length, ier)

  end subroutine

  subroutine iricmi_rin_integer(name, val, ier)
    character(*), intent(in):: name
    integer, intent(in):: val
    integer, intent(out):: ier

    call iricmi_rin_integer_f2c &
      (name, val, ier)

  end subroutine

  subroutine iricmi_rin_real(name, val, ier)
    character(*), intent(in):: name
    double precision, intent(in):: val
    integer, intent(out):: ier

    call iricmi_rin_real_f2c &
      (name, val, ier)

  end subroutine

  subroutine iricmi_rout_integer(name, val, ier)
    character(*), intent(in):: name
    integer, intent(in):: val
    integer, intent(out):: ier

    call iricmi_rout_integer_f2c &
      (name, val, ier)

  end subroutine

  subroutine iricmi_rout_real(name, val, ier)
    character(*), intent(in):: name
    double precision, intent(in):: val
    integer, intent(out):: ier

    call iricmi_rout_real_f2c &
      (name, val, ier)

  end subroutine



  ! from iricmi_complex.h

  subroutine iricmi_read_grid_complex_node_withgridid_1d(gid, groupname, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: groupname
    integer, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_complex_node_withgridid_f2c &
      (gid, groupname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_complex_node_withgridid_2d(gid, groupname, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: groupname
    integer, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_complex_node_withgridid_f2c &
      (gid, groupname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_complex_node_withgridid_3d(gid, groupname, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: groupname
    integer, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_complex_node_withgridid_f2c &
      (gid, groupname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_complex_cell_withgridid_1d(gid, groupname, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: groupname
    integer, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_complex_cell_withgridid_f2c &
      (gid, groupname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_complex_cell_withgridid_2d(gid, groupname, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: groupname
    integer, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_complex_cell_withgridid_f2c &
      (gid, groupname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_complex_cell_withgridid_3d(gid, groupname, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: groupname
    integer, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_complex_cell_withgridid_f2c &
      (gid, groupname, v_arr, ier)

  end subroutine



  ! from iricmi_geo.h

  subroutine iricmi_geo_riversurvey_open(filename, id, ier)
    character(*), intent(in):: filename
    integer, intent(out):: id
    integer, intent(out):: ier

    call iricmi_geo_riversurvey_open_f2c &
      (filename, id, ier)

  end subroutine

  subroutine iricmi_geo_riversurvey_read_count(id, count, ier)
    integer, intent(in):: id
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_geo_riversurvey_read_count_f2c &
      (id, count, ier)

  end subroutine

  subroutine iricmi_geo_riversurvey_read_position(id, pointid, x, y, ier)
    integer, intent(in):: id
    integer, intent(in):: pointid
    double precision, intent(out):: x
    double precision, intent(out):: y
    integer, intent(out):: ier

    call iricmi_geo_riversurvey_read_position_f2c &
      (id, pointid, x, y, ier)

  end subroutine

  subroutine iricmi_geo_riversurvey_read_direction(id, pointid, vx, vy, ier)
    integer, intent(in):: id
    integer, intent(in):: pointid
    double precision, intent(out):: vx
    double precision, intent(out):: vy
    integer, intent(out):: ier

    call iricmi_geo_riversurvey_read_direction_f2c &
      (id, pointid, vx, vy, ier)

  end subroutine

  subroutine iricmi_geo_riversurvey_read_name(id, pointid, strvalue, ier)
    integer, intent(in):: id
    integer, intent(in):: pointid
    character(*), intent(out):: strvalue
    integer, intent(out):: ier

    call iricmi_geo_riversurvey_read_name_f2c &
      (id, pointid, strvalue, ier)

  end subroutine

  subroutine iricmi_geo_riversurvey_read_realname(id, pointid, name, ier)
    integer, intent(in):: id
    integer, intent(in):: pointid
    double precision, intent(out):: name
    integer, intent(out):: ier

    call iricmi_geo_riversurvey_read_realname_f2c &
      (id, pointid, name, ier)

  end subroutine

  subroutine iricmi_geo_riversurvey_read_leftshift(id, pointid, shift, ier)
    integer, intent(in):: id
    integer, intent(in):: pointid
    double precision, intent(out):: shift
    integer, intent(out):: ier

    call iricmi_geo_riversurvey_read_leftshift_f2c &
      (id, pointid, shift, ier)

  end subroutine

  subroutine iricmi_geo_riversurvey_read_altitudecount(id, pointid, count, ier)
    integer, intent(in):: id
    integer, intent(in):: pointid
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_geo_riversurvey_read_altitudecount_f2c &
      (id, pointid, count, ier)

  end subroutine

  subroutine iricmi_geo_riversurvey_read_altitudes(id, pointid, position_arr, height_arr, active_arr, ier)
    integer, intent(in):: id
    integer, intent(in):: pointid
    double precision, dimension(:), intent(out):: position_arr
    double precision, dimension(:), intent(out):: height_arr
    integer, dimension(:), intent(out):: active_arr
    integer, intent(out):: ier

    call iricmi_geo_riversurvey_read_altitudes_f2c &
      (id, pointid, position_arr, height_arr, active_arr, ier)

  end subroutine

  subroutine iricmi_geo_riversurvey_read_fixedpointl(id, pointid, set, directionX, directionY, index, ier)
    integer, intent(in):: id
    integer, intent(in):: pointid
    integer, intent(out):: set
    double precision, intent(out):: directionX
    double precision, intent(out):: directionY
    integer, intent(out):: index
    integer, intent(out):: ier

    call iricmi_geo_riversurvey_read_fixedpointl_f2c &
      (id, pointid, set, directionX, directionY, index, ier)

  end subroutine

  subroutine iricmi_geo_riversurvey_read_fixedpointr(id, pointid, set, directionX, directionY, index, ier)
    integer, intent(in):: id
    integer, intent(in):: pointid
    integer, intent(out):: set
    double precision, intent(out):: directionX
    double precision, intent(out):: directionY
    integer, intent(out):: index
    integer, intent(out):: ier

    call iricmi_geo_riversurvey_read_fixedpointr_f2c &
      (id, pointid, set, directionX, directionY, index, ier)

  end subroutine

  subroutine iricmi_geo_riversurvey_read_watersurfaceelevation(id, pointid, set, value, ier)
    integer, intent(in):: id
    integer, intent(in):: pointid
    integer, intent(out):: set
    double precision, intent(out):: value
    integer, intent(out):: ier

    call iricmi_geo_riversurvey_read_watersurfaceelevation_f2c &
      (id, pointid, set, value, ier)

  end subroutine

  subroutine iricmi_geo_riversurvey_close(id, ier)
    integer, intent(in):: id
    integer, intent(out):: ier

    call iricmi_geo_riversurvey_close_f2c &
      (id, ier)

  end subroutine



  ! from iricmi_geoutil.h

  subroutine iricmi_read_geo_count(name, count, ier)
    character(*), intent(in):: name
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_read_geo_count_f2c &
      (name, count, ier)

  end subroutine

  subroutine iricmi_read_geo_filename(name, geoid, strvalue, type, ier)
    character(*), intent(in):: name
    integer, intent(in):: geoid
    character(*), intent(out):: strvalue
    integer, intent(out):: type
    integer, intent(out):: ier

    call iricmi_read_geo_filename_f2c &
      (name, geoid, strvalue, type, ier)

  end subroutine



  ! from iricmi_grid.h

  subroutine iricmi_read_grid2d_str_size_withgridid(gid, isize, jsize, ier)
    integer, intent(in):: gid
    integer, intent(out):: isize
    integer, intent(out):: jsize
    integer, intent(out):: ier

    call iricmi_read_grid2d_str_size_withgridid_f2c &
      (gid, isize, jsize, ier)

  end subroutine

  subroutine iricmi_read_grid2d_coords_withgridid_1d(gid, x_arr, y_arr, ier)
    integer, intent(in):: gid
    double precision, dimension(:), intent(out):: x_arr
    double precision, dimension(:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_grid2d_coords_withgridid_f2c &
      (gid, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_grid2d_coords_withgridid_2d(gid, x_arr, y_arr, ier)
    integer, intent(in):: gid
    double precision, dimension(:,:), intent(out):: x_arr
    double precision, dimension(:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_grid2d_coords_withgridid_f2c &
      (gid, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_grid3d_str_size_withgridid(gid, isize, jsize, ksize, ier)
    integer, intent(in):: gid
    integer, intent(out):: isize
    integer, intent(out):: jsize
    integer, intent(out):: ksize
    integer, intent(out):: ier

    call iricmi_read_grid3d_str_size_withgridid_f2c &
      (gid, isize, jsize, ksize, ier)

  end subroutine

  subroutine iricmi_read_grid3d_coords_withgridid_1d(gid, x_arr, y_arr, z_arr, ier)
    integer, intent(in):: gid
    double precision, dimension(:), intent(out):: x_arr
    double precision, dimension(:), intent(out):: y_arr
    double precision, dimension(:), intent(out):: z_arr
    integer, intent(out):: ier

    call iricmi_read_grid3d_coords_withgridid_f2c &
      (gid, x_arr, y_arr, z_arr, ier)

  end subroutine

  subroutine iricmi_read_grid3d_coords_withgridid_3d(gid, x_arr, y_arr, z_arr, ier)
    integer, intent(in):: gid
    double precision, dimension(:,:,:), intent(out):: x_arr
    double precision, dimension(:,:,:), intent(out):: y_arr
    double precision, dimension(:,:,:), intent(out):: z_arr
    integer, intent(out):: ier

    call iricmi_read_grid3d_coords_withgridid_f2c &
      (gid, x_arr, y_arr, z_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_triangleelementssize_withgridid(gid, size, ier)
    integer, intent(in):: gid
    integer, intent(out):: size
    integer, intent(out):: ier

    call iricmi_read_grid_triangleelementssize_withgridid_f2c &
      (gid, size, ier)

  end subroutine

  subroutine iricmi_read_grid_triangleelements_withgridid(gid, id_arr, ier)
    integer, intent(in):: gid
    integer, dimension(:), intent(out):: id_arr
    integer, intent(out):: ier

    call iricmi_read_grid_triangleelements_withgridid_f2c &
      (gid, id_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_nodecount_withgridid(gid, count, ier)
    integer, intent(in):: gid
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_read_grid_nodecount_withgridid_f2c &
      (gid, count, ier)

  end subroutine

  subroutine iricmi_read_grid_cellcount_withgridid(gid, count, ier)
    integer, intent(in):: gid
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_read_grid_cellcount_withgridid_f2c &
      (gid, count, ier)

  end subroutine

  subroutine iricmi_read_grid_ifacecount_withgridid(gid, count, ier)
    integer, intent(in):: gid
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_read_grid_ifacecount_withgridid_f2c &
      (gid, count, ier)

  end subroutine

  subroutine iricmi_read_grid_jfacecount_withgridid(gid, count, ier)
    integer, intent(in):: gid
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_read_grid_jfacecount_withgridid_f2c &
      (gid, count, ier)

  end subroutine

  subroutine iricmi_read_grid_kfacecount_withgridid(gid, count, ier)
    integer, intent(in):: gid
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_read_grid_kfacecount_withgridid_f2c &
      (gid, count, ier)

  end subroutine

  subroutine iricmi_read_grid_real_node_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_real_node_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_real_node_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_real_node_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_real_node_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_real_node_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_integer_node_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_integer_node_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_integer_node_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_integer_node_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_integer_node_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_integer_node_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_real_cell_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_real_cell_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_real_cell_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_real_cell_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_real_cell_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_real_cell_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_integer_cell_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_integer_cell_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_integer_cell_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_integer_cell_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_integer_cell_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_integer_cell_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functionaldimensionsize_withgridid(gid, name, dimname, count, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    character(*), intent(in):: dimname
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_read_grid_functionaldimensionsize_withgridid_f2c &
      (gid, name, dimname, count, ier)

  end subroutine

  subroutine iricmi_read_grid_functionaldimension_integer_withgridid(gid, name, dimname, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    character(*), intent(in):: dimname
    integer, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functionaldimension_integer_withgridid_f2c &
      (gid, name, dimname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functionaldimension_real_withgridid(gid, name, dimname, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    character(*), intent(in):: dimname
    double precision, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functionaldimension_real_withgridid_f2c &
      (gid, name, dimname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functionaltimesize_withgridid(gid, name, count, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_read_grid_functionaltimesize_withgridid_f2c &
      (gid, name, count, ier)

  end subroutine

  subroutine iricmi_read_grid_functionaltime_withgridid(gid, name, time_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:), intent(out):: time_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functionaltime_withgridid_f2c &
      (gid, name, time_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_integer_node_withgridid_1d(gid, name, dimid, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, intent(in):: dimid
    integer, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_integer_node_withgridid_f2c &
      (gid, name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_integer_node_withgridid_2d(gid, name, dimid, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, intent(in):: dimid
    integer, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_integer_node_withgridid_f2c &
      (gid, name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_integer_node_withgridid_3d(gid, name, dimid, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, intent(in):: dimid
    integer, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_integer_node_withgridid_f2c &
      (gid, name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_real_node_withgridid_1d(gid, name, dimid, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, intent(in):: dimid
    double precision, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_real_node_withgridid_f2c &
      (gid, name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_real_node_withgridid_2d(gid, name, dimid, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, intent(in):: dimid
    double precision, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_real_node_withgridid_f2c &
      (gid, name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_real_node_withgridid_3d(gid, name, dimid, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, intent(in):: dimid
    double precision, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_real_node_withgridid_f2c &
      (gid, name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_integer_cell_withgridid_1d(gid, name, dimid, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, intent(in):: dimid
    integer, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_integer_cell_withgridid_f2c &
      (gid, name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_integer_cell_withgridid_2d(gid, name, dimid, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, intent(in):: dimid
    integer, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_integer_cell_withgridid_f2c &
      (gid, name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_integer_cell_withgridid_3d(gid, name, dimid, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, intent(in):: dimid
    integer, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_integer_cell_withgridid_f2c &
      (gid, name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_real_cell_withgridid_1d(gid, name, dimid, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, intent(in):: dimid
    double precision, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_real_cell_withgridid_f2c &
      (gid, name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_real_cell_withgridid_2d(gid, name, dimid, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, intent(in):: dimid
    double precision, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_real_cell_withgridid_f2c &
      (gid, name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_real_cell_withgridid_3d(gid, name, dimid, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, intent(in):: dimid
    double precision, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_real_cell_withgridid_f2c &
      (gid, name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_write_grid1d_coords_withgridid(isize, x_arr, gid, ier)
    integer, intent(in):: isize
    double precision, dimension(:), intent(in):: x_arr
    integer, intent(out):: gid
    integer, intent(out):: ier

    call iricmi_write_grid1d_coords_withgridid_f2c &
      (isize, x_arr, gid, ier)

  end subroutine

  subroutine iricmi_write_grid2d_coords_withgridid_1d(isize, jsize, x_arr, y_arr, gid, ier)
    integer, intent(in):: isize
    integer, intent(in):: jsize
    double precision, dimension(:), intent(in):: x_arr
    double precision, dimension(:), intent(in):: y_arr
    integer, intent(out):: gid
    integer, intent(out):: ier

    call iricmi_write_grid2d_coords_withgridid_f2c &
      (isize, jsize, x_arr, y_arr, gid, ier)

  end subroutine

  subroutine iricmi_write_grid2d_coords_withgridid_2d(isize, jsize, x_arr, y_arr, gid, ier)
    integer, intent(in):: isize
    integer, intent(in):: jsize
    double precision, dimension(:,:), intent(in):: x_arr
    double precision, dimension(:,:), intent(in):: y_arr
    integer, intent(out):: gid
    integer, intent(out):: ier

    call iricmi_write_grid2d_coords_withgridid_f2c &
      (isize, jsize, x_arr, y_arr, gid, ier)

  end subroutine

  subroutine iricmi_write_grid3d_coords_withgridid_1d(isize, jsize, ksize, x_arr, y_arr, z_arr, gid, ier)
    integer, intent(in):: isize
    integer, intent(in):: jsize
    integer, intent(in):: ksize
    double precision, dimension(:), intent(in):: x_arr
    double precision, dimension(:), intent(in):: y_arr
    double precision, dimension(:), intent(in):: z_arr
    integer, intent(out):: gid
    integer, intent(out):: ier

    call iricmi_write_grid3d_coords_withgridid_f2c &
      (isize, jsize, ksize, x_arr, y_arr, z_arr, gid, ier)

  end subroutine

  subroutine iricmi_write_grid3d_coords_withgridid_3d(isize, jsize, ksize, x_arr, y_arr, z_arr, gid, ier)
    integer, intent(in):: isize
    integer, intent(in):: jsize
    integer, intent(in):: ksize
    double precision, dimension(:,:,:), intent(in):: x_arr
    double precision, dimension(:,:,:), intent(in):: y_arr
    double precision, dimension(:,:,:), intent(in):: z_arr
    integer, intent(out):: gid
    integer, intent(out):: ier

    call iricmi_write_grid3d_coords_withgridid_f2c &
      (isize, jsize, ksize, x_arr, y_arr, z_arr, gid, ier)

  end subroutine

  subroutine iricmi_write_namedgrid1d_coords_withgridid(name, isize, x_arr, gid, ier)
    character(*), intent(in):: name
    integer, intent(in):: isize
    double precision, dimension(:), intent(in):: x_arr
    integer, intent(out):: gid
    integer, intent(out):: ier

    call iricmi_write_namedgrid1d_coords_withgridid_f2c &
      (name, isize, x_arr, gid, ier)

  end subroutine

  subroutine iricmi_write_namedgrid2d_coords_withgridid_1d(name, isize, jsize, x_arr, y_arr, gid, ier)
    character(*), intent(in):: name
    integer, intent(in):: isize
    integer, intent(in):: jsize
    double precision, dimension(:), intent(in):: x_arr
    double precision, dimension(:), intent(in):: y_arr
    integer, intent(out):: gid
    integer, intent(out):: ier

    call iricmi_write_namedgrid2d_coords_withgridid_f2c &
      (name, isize, jsize, x_arr, y_arr, gid, ier)

  end subroutine

  subroutine iricmi_write_namedgrid2d_coords_withgridid_2d(name, isize, jsize, x_arr, y_arr, gid, ier)
    character(*), intent(in):: name
    integer, intent(in):: isize
    integer, intent(in):: jsize
    double precision, dimension(:,:), intent(in):: x_arr
    double precision, dimension(:,:), intent(in):: y_arr
    integer, intent(out):: gid
    integer, intent(out):: ier

    call iricmi_write_namedgrid2d_coords_withgridid_f2c &
      (name, isize, jsize, x_arr, y_arr, gid, ier)

  end subroutine

  subroutine iricmi_write_namedgrid3d_coords_withgridid_1d(name, isize, jsize, ksize, x_arr, y_arr, z_arr, gid, ier)
    character(*), intent(in):: name
    integer, intent(in):: isize
    integer, intent(in):: jsize
    integer, intent(in):: ksize
    double precision, dimension(:), intent(in):: x_arr
    double precision, dimension(:), intent(in):: y_arr
    double precision, dimension(:), intent(in):: z_arr
    integer, intent(out):: gid
    integer, intent(out):: ier

    call iricmi_write_namedgrid3d_coords_withgridid_f2c &
      (name, isize, jsize, ksize, x_arr, y_arr, z_arr, gid, ier)

  end subroutine

  subroutine iricmi_write_namedgrid3d_coords_withgridid_3d(name, isize, jsize, ksize, x_arr, y_arr, z_arr, gid, ier)
    character(*), intent(in):: name
    integer, intent(in):: isize
    integer, intent(in):: jsize
    integer, intent(in):: ksize
    double precision, dimension(:,:,:), intent(in):: x_arr
    double precision, dimension(:,:,:), intent(in):: y_arr
    double precision, dimension(:,:,:), intent(in):: z_arr
    integer, intent(out):: gid
    integer, intent(out):: ier

    call iricmi_write_namedgrid3d_coords_withgridid_f2c &
      (name, isize, jsize, ksize, x_arr, y_arr, z_arr, gid, ier)

  end subroutine



  ! from iricmi_grid_solverlib.h

  subroutine iricmi_read_grid2d_open_withgridid(gid, grid_handle, ier)
    integer, intent(in):: gid
    integer, intent(out):: grid_handle
    integer, intent(out):: ier

    call iricmi_read_grid2d_open_withgridid_f2c &
      (gid, grid_handle, ier)

  end subroutine

  subroutine iricmi_read_grid2d_close(grid_handle, ier)
    integer, intent(in):: grid_handle
    integer, intent(out):: ier

    call iricmi_read_grid2d_close_f2c &
      (grid_handle, ier)

  end subroutine

  subroutine iricmi_read_grid2d_cellarea(grid_handle, cellId, area, ier)
    integer, intent(in):: grid_handle
    integer, intent(in):: cellId
    double precision, intent(out):: area
    integer, intent(out):: ier

    call iricmi_read_grid2d_cellarea_f2c &
      (grid_handle, cellId, area, ier)

  end subroutine

  subroutine iricmi_read_grid2d_findcell(grid_handle, x, y, cellId, ier)
    integer, intent(in):: grid_handle
    double precision, intent(in):: x
    double precision, intent(in):: y
    integer, intent(out):: cellId
    integer, intent(out):: ier

    call iricmi_read_grid2d_findcell_f2c &
      (grid_handle, x, y, cellId, ier)

  end subroutine

  subroutine iricmi_read_grid2d_cellnodecount(grid_handle, cellId, count, ier)
    integer, intent(in):: grid_handle
    integer, intent(in):: cellId
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_read_grid2d_cellnodecount_f2c &
      (grid_handle, cellId, count, ier)

  end subroutine

  subroutine iricmi_read_grid2d_interpolate(grid_handle, x, y, ok, count, nodeids_arr, weights_arr, ier)
    integer, intent(in):: grid_handle
    double precision, intent(in):: x
    double precision, intent(in):: y
    integer, intent(out):: ok
    integer, intent(out):: count
    integer, dimension(:), intent(out):: nodeids_arr
    double precision, dimension(:), intent(out):: weights_arr
    integer, intent(out):: ier

    call iricmi_read_grid2d_interpolate_f2c &
      (grid_handle, x, y, ok, count, nodeids_arr, weights_arr, ier)

  end subroutine

  subroutine iricmi_read_grid2d_interpolatewithcell(grid_handle, x, y, cellId, nodeids_arr, weights_arr, ier)
    integer, intent(in):: grid_handle
    double precision, intent(in):: x
    double precision, intent(in):: y
    integer, intent(in):: cellId
    integer, dimension(:), intent(out):: nodeids_arr
    double precision, dimension(:), intent(out):: weights_arr
    integer, intent(out):: ier

    call iricmi_read_grid2d_interpolatewithcell_f2c &
      (grid_handle, x, y, cellId, nodeids_arr, weights_arr, ier)

  end subroutine



  ! from iricmi_gui_coorp.h

  subroutine iricmi_check_cancel(ier)
    integer, intent(out):: ier

    call iricmi_check_cancel_f2c &
      (ier)

  end subroutine

  subroutine iricmi_check_update(ier)
    integer, intent(out):: ier

    call iricmi_check_update_f2c &
      (ier)

  end subroutine



  ! from iricmi_model.h

  subroutine iricmi_model_init(ier)
    integer, intent(out):: ier

    call iricmi_model_init_f2c &
      (ier)

  end subroutine

  subroutine iricmi_model_terminate(ier)
    integer, intent(out):: ier

    call iricmi_model_terminate_f2c &
      (ier)

  end subroutine

  subroutine iricmi_model_sync(ier)
    integer, intent(out):: ier

    call iricmi_model_sync_f2c &
      (ier)

  end subroutine

  subroutine iricmi_model_dump(ier)
    integer, intent(out):: ier

    call iricmi_model_dump_f2c &
      (ier)

  end subroutine



  ! from iricmi_not_withbaseid.h

  subroutine iricmi_read_complex_count(groupname, num, ier)
    character(*), intent(in):: groupname
    integer, intent(out):: num
    integer, intent(out):: ier

    call iricmi_read_complex_count_f2c &
      (groupname, num, ier)

  end subroutine

  subroutine iricmi_read_complex_integer(groupname, num, name, value, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    integer, intent(out):: value
    integer, intent(out):: ier

    call iricmi_read_complex_integer_f2c &
      (groupname, num, name, value, ier)

  end subroutine

  subroutine iricmi_read_complex_real(groupname, num, name, value, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, intent(out):: value
    integer, intent(out):: ier

    call iricmi_read_complex_real_f2c &
      (groupname, num, name, value, ier)

  end subroutine

  subroutine iricmi_read_complex_realsingle(groupname, num, name, value, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    real, intent(out):: value
    integer, intent(out):: ier

    call iricmi_read_complex_realsingle_f2c &
      (groupname, num, name, value, ier)

  end subroutine

  subroutine iricmi_read_complex_stringlen(groupname, num, name, length, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    integer, intent(out):: length
    integer, intent(out):: ier

    call iricmi_read_complex_stringlen_f2c &
      (groupname, num, name, length, ier)

  end subroutine

  subroutine iricmi_read_complex_string(groupname, num, name, strvalue, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(out):: strvalue
    integer, intent(out):: ier

    call iricmi_read_complex_string_f2c &
      (groupname, num, name, strvalue, ier)

  end subroutine

  subroutine iricmi_read_complex_functionalsize(groupname, num, name, size, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    integer, intent(out):: size
    integer, intent(out):: ier

    call iricmi_read_complex_functionalsize_f2c &
      (groupname, num, name, size, ier)

  end subroutine

  subroutine iricmi_read_complex_functional_1d(groupname, num, name, x_arr, y_arr, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, dimension(:), intent(out):: x_arr
    double precision, dimension(:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_complex_functional_f2c &
      (groupname, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_complex_functional_2d(groupname, num, name, x_arr, y_arr, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(out):: x_arr
    double precision, dimension(:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_complex_functional_f2c &
      (groupname, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_complex_functional_3d(groupname, num, name, x_arr, y_arr, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(out):: x_arr
    double precision, dimension(:,:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_complex_functional_f2c &
      (groupname, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_complex_functional_4d(groupname, num, name, x_arr, y_arr, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, dimension(:,:,:,:), intent(out):: x_arr
    double precision, dimension(:,:,:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_complex_functional_f2c &
      (groupname, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_complex_functionalwithname_1d(groupname, num, name, paramname, v_arr, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    double precision, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_complex_functionalwithname_f2c &
      (groupname, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_complex_functionalwithname_2d(groupname, num, name, paramname, v_arr, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    double precision, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_complex_functionalwithname_f2c &
      (groupname, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_complex_functionalwithname_3d(groupname, num, name, paramname, v_arr, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    double precision, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_complex_functionalwithname_f2c &
      (groupname, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_complex_functionalwithname_4d(groupname, num, name, paramname, v_arr, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    double precision, dimension(:,:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_complex_functionalwithname_f2c &
      (groupname, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_complex_functional_realsingle_1d(groupname, num, name, x_arr, y_arr, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    real, dimension(:), intent(out):: x_arr
    real, dimension(:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_complex_functional_realsingle_f2c &
      (groupname, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_complex_functional_realsingle_2d(groupname, num, name, x_arr, y_arr, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    real, dimension(:,:), intent(out):: x_arr
    real, dimension(:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_complex_functional_realsingle_f2c &
      (groupname, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_complex_functional_realsingle_3d(groupname, num, name, x_arr, y_arr, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    real, dimension(:,:,:), intent(out):: x_arr
    real, dimension(:,:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_complex_functional_realsingle_f2c &
      (groupname, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_complex_functional_realsingle_4d(groupname, num, name, x_arr, y_arr, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    real, dimension(:,:,:,:), intent(out):: x_arr
    real, dimension(:,:,:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_complex_functional_realsingle_f2c &
      (groupname, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_complex_functionalwithname_realsingle_1d(groupname, num, name, paramname, v_arr, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    real, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_complex_functionalwithname_realsingle_f2c &
      (groupname, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_complex_functionalwithname_realsingle_2d(groupname, num, name, paramname, v_arr, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    real, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_complex_functionalwithname_realsingle_f2c &
      (groupname, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_complex_functionalwithname_realsingle_3d(groupname, num, name, paramname, v_arr, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    real, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_complex_functionalwithname_realsingle_f2c &
      (groupname, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_complex_functionalwithname_realsingle_4d(groupname, num, name, paramname, v_arr, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    real, dimension(:,:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_complex_functionalwithname_realsingle_f2c &
      (groupname, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_complex_functionalwithname_stringlen(groupname, num, name, paramname, length, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    integer, intent(out):: length
    integer, intent(out):: ier

    call iricmi_read_complex_functionalwithname_stringlen_f2c &
      (groupname, num, name, paramname, length, ier)

  end subroutine

  subroutine iricmi_read_complex_functionalwithname_string(groupname, num, name, paramname, strvalue, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    character(*), intent(out):: strvalue
    integer, intent(out):: ier

    call iricmi_read_complex_functionalwithname_string_f2c &
      (groupname, num, name, paramname, strvalue, ier)

  end subroutine

  subroutine iricmi_rin_complex_integer(groupname, num, name, val, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    integer, intent(in):: val
    integer, intent(out):: ier

    call iricmi_rin_complex_integer_f2c &
      (groupname, num, name, val, ier)

  end subroutine

  subroutine iricmi_rin_complex_real(groupname, num, name, val, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, intent(in):: val
    integer, intent(out):: ier

    call iricmi_rin_complex_real_f2c &
      (groupname, num, name, val, ier)

  end subroutine

  subroutine iricmi_rout_complex_integer(groupname, num, name, val, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    integer, intent(in):: val
    integer, intent(out):: ier

    call iricmi_rout_complex_integer_f2c &
      (groupname, num, name, val, ier)

  end subroutine

  subroutine iricmi_rout_complex_real(groupname, num, name, val, ier)
    character(*), intent(in):: groupname
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, intent(in):: val
    integer, intent(out):: ier

    call iricmi_rout_complex_real_f2c &
      (groupname, num, name, val, ier)

  end subroutine



  ! from iricmi_not_withgridid.h

  subroutine iricmi_read_bc_count(type, num, ier)
    character(*), intent(in):: type
    integer, intent(out):: num
    integer, intent(out):: ier

    call iricmi_read_bc_count_f2c &
      (type, num, ier)

  end subroutine

  subroutine iricmi_read_bc_indicessize(type, num, size, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    integer, intent(out):: size
    integer, intent(out):: ier

    call iricmi_read_bc_indicessize_f2c &
      (type, num, size, ier)

  end subroutine

  subroutine iricmi_read_bc_indices_1d(type, num, idx_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    integer, dimension(:), intent(out):: idx_arr
    integer, intent(out):: ier

    call iricmi_read_bc_indices_f2c &
      (type, num, idx_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_indices_2d(type, num, idx_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    integer, dimension(:,:), intent(out):: idx_arr
    integer, intent(out):: ier

    call iricmi_read_bc_indices_f2c &
      (type, num, idx_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_indices_3d(type, num, idx_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    integer, dimension(:,:,:), intent(out):: idx_arr
    integer, intent(out):: ier

    call iricmi_read_bc_indices_f2c &
      (type, num, idx_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_integer(type, num, name, value, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    integer, intent(out):: value
    integer, intent(out):: ier

    call iricmi_read_bc_integer_f2c &
      (type, num, name, value, ier)

  end subroutine

  subroutine iricmi_read_bc_real(type, num, name, value, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, intent(out):: value
    integer, intent(out):: ier

    call iricmi_read_bc_real_f2c &
      (type, num, name, value, ier)

  end subroutine

  subroutine iricmi_read_bc_realsingle(type, num, name, value, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    real, intent(out):: value
    integer, intent(out):: ier

    call iricmi_read_bc_realsingle_f2c &
      (type, num, name, value, ier)

  end subroutine

  subroutine iricmi_read_bc_stringlen(type, num, name, length, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    integer, intent(out):: length
    integer, intent(out):: ier

    call iricmi_read_bc_stringlen_f2c &
      (type, num, name, length, ier)

  end subroutine

  subroutine iricmi_read_bc_string(type, num, name, strvalue, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(out):: strvalue
    integer, intent(out):: ier

    call iricmi_read_bc_string_f2c &
      (type, num, name, strvalue, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalsize(type, num, name, size, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    integer, intent(out):: size
    integer, intent(out):: ier

    call iricmi_read_bc_functionalsize_f2c &
      (type, num, name, size, ier)

  end subroutine

  subroutine iricmi_read_bc_functional_1d(type, num, name, x_arr, y_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, dimension(:), intent(out):: x_arr
    double precision, dimension(:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functional_f2c &
      (type, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functional_2d(type, num, name, x_arr, y_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(out):: x_arr
    double precision, dimension(:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functional_f2c &
      (type, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functional_3d(type, num, name, x_arr, y_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(out):: x_arr
    double precision, dimension(:,:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functional_f2c &
      (type, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functional_4d(type, num, name, x_arr, y_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, dimension(:,:,:,:), intent(out):: x_arr
    double precision, dimension(:,:,:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functional_f2c &
      (type, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_1d(type, num, name, paramname, v_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    double precision, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_f2c &
      (type, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_2d(type, num, name, paramname, v_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    double precision, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_f2c &
      (type, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_3d(type, num, name, paramname, v_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    double precision, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_f2c &
      (type, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_4d(type, num, name, paramname, v_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    double precision, dimension(:,:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_f2c &
      (type, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functional_realsingle_1d(type, num, name, x_arr, y_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    real, dimension(:), intent(out):: x_arr
    real, dimension(:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functional_realsingle_f2c &
      (type, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functional_realsingle_2d(type, num, name, x_arr, y_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    real, dimension(:,:), intent(out):: x_arr
    real, dimension(:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functional_realsingle_f2c &
      (type, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functional_realsingle_3d(type, num, name, x_arr, y_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    real, dimension(:,:,:), intent(out):: x_arr
    real, dimension(:,:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functional_realsingle_f2c &
      (type, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functional_realsingle_4d(type, num, name, x_arr, y_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    real, dimension(:,:,:,:), intent(out):: x_arr
    real, dimension(:,:,:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functional_realsingle_f2c &
      (type, num, name, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_realsingle_1d(type, num, name, paramname, v_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    real, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_realsingle_f2c &
      (type, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_realsingle_2d(type, num, name, paramname, v_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    real, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_realsingle_f2c &
      (type, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_realsingle_3d(type, num, name, paramname, v_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    real, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_realsingle_f2c &
      (type, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_realsingle_4d(type, num, name, paramname, v_arr, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    real, dimension(:,:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_realsingle_f2c &
      (type, num, name, paramname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_stringlen(type, num, name, paramname, length, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    integer, intent(out):: length
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_stringlen_f2c &
      (type, num, name, paramname, length, ier)

  end subroutine

  subroutine iricmi_read_bc_functionalwithname_string(type, num, name, paramname, strvalue, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    character(*), intent(in):: paramname
    character(*), intent(out):: strvalue
    integer, intent(out):: ier

    call iricmi_read_bc_functionalwithname_string_f2c &
      (type, num, name, paramname, strvalue, ier)

  end subroutine

  subroutine iricmi_rin_bc_integer(type, num, name, val, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    integer, intent(in):: val
    integer, intent(out):: ier

    call iricmi_rin_bc_integer_f2c &
      (type, num, name, val, ier)

  end subroutine

  subroutine iricmi_rin_bc_real(type, num, name, val, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, intent(in):: val
    integer, intent(out):: ier

    call iricmi_rin_bc_real_f2c &
      (type, num, name, val, ier)

  end subroutine

  subroutine iricmi_rout_bc_integer(type, num, name, val, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    integer, intent(in):: val
    integer, intent(out):: ier

    call iricmi_rout_bc_integer_f2c &
      (type, num, name, val, ier)

  end subroutine

  subroutine iricmi_rout_bc_real(type, num, name, val, ier)
    character(*), intent(in):: type
    integer, intent(in):: num
    character(*), intent(in):: name
    double precision, intent(in):: val
    integer, intent(out):: ier

    call iricmi_rout_bc_real_f2c &
      (type, num, name, val, ier)

  end subroutine

  subroutine iricmi_read_grid_complex_node_1d(groupname, v_arr, ier)
    character(*), intent(in):: groupname
    integer, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_complex_node_f2c &
      (groupname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_complex_node_2d(groupname, v_arr, ier)
    character(*), intent(in):: groupname
    integer, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_complex_node_f2c &
      (groupname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_complex_node_3d(groupname, v_arr, ier)
    character(*), intent(in):: groupname
    integer, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_complex_node_f2c &
      (groupname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_complex_cell_1d(groupname, v_arr, ier)
    character(*), intent(in):: groupname
    integer, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_complex_cell_f2c &
      (groupname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_complex_cell_2d(groupname, v_arr, ier)
    character(*), intent(in):: groupname
    integer, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_complex_cell_f2c &
      (groupname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_complex_cell_3d(groupname, v_arr, ier)
    character(*), intent(in):: groupname
    integer, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_complex_cell_f2c &
      (groupname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid2d_str_size(isize, jsize, ier)
    integer, intent(out):: isize
    integer, intent(out):: jsize
    integer, intent(out):: ier

    call iricmi_read_grid2d_str_size_f2c &
      (isize, jsize, ier)

  end subroutine

  subroutine iricmi_read_grid2d_coords_1d(x_arr, y_arr, ier)
    double precision, dimension(:), intent(out):: x_arr
    double precision, dimension(:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_grid2d_coords_f2c &
      (x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_grid2d_coords_2d(x_arr, y_arr, ier)
    double precision, dimension(:,:), intent(out):: x_arr
    double precision, dimension(:,:), intent(out):: y_arr
    integer, intent(out):: ier

    call iricmi_read_grid2d_coords_f2c &
      (x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_read_grid3d_str_size(isize, jsize, ksize, ier)
    integer, intent(out):: isize
    integer, intent(out):: jsize
    integer, intent(out):: ksize
    integer, intent(out):: ier

    call iricmi_read_grid3d_str_size_f2c &
      (isize, jsize, ksize, ier)

  end subroutine

  subroutine iricmi_read_grid3d_coords_1d(x_arr, y_arr, z_arr, ier)
    double precision, dimension(:), intent(out):: x_arr
    double precision, dimension(:), intent(out):: y_arr
    double precision, dimension(:), intent(out):: z_arr
    integer, intent(out):: ier

    call iricmi_read_grid3d_coords_f2c &
      (x_arr, y_arr, z_arr, ier)

  end subroutine

  subroutine iricmi_read_grid3d_coords_3d(x_arr, y_arr, z_arr, ier)
    double precision, dimension(:,:,:), intent(out):: x_arr
    double precision, dimension(:,:,:), intent(out):: y_arr
    double precision, dimension(:,:,:), intent(out):: z_arr
    integer, intent(out):: ier

    call iricmi_read_grid3d_coords_f2c &
      (x_arr, y_arr, z_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_triangleelementssize(size, ier)
    integer, intent(out):: size
    integer, intent(out):: ier

    call iricmi_read_grid_triangleelementssize_f2c &
      (size, ier)

  end subroutine

  subroutine iricmi_read_grid_triangleelements(id_arr, ier)
    integer, dimension(:), intent(out):: id_arr
    integer, intent(out):: ier

    call iricmi_read_grid_triangleelements_f2c &
      (id_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_nodecount(count, ier)
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_read_grid_nodecount_f2c &
      (count, ier)

  end subroutine

  subroutine iricmi_read_grid_cellcount(count, ier)
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_read_grid_cellcount_f2c &
      (count, ier)

  end subroutine

  subroutine iricmi_read_grid_ifacecount(count, ier)
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_read_grid_ifacecount_f2c &
      (count, ier)

  end subroutine

  subroutine iricmi_read_grid_jfacecount(count, ier)
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_read_grid_jfacecount_f2c &
      (count, ier)

  end subroutine

  subroutine iricmi_read_grid_kfacecount(count, ier)
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_read_grid_kfacecount_f2c &
      (count, ier)

  end subroutine

  subroutine iricmi_read_grid_real_node_1d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_real_node_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_real_node_2d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_real_node_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_real_node_3d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_real_node_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_integer_node_1d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_integer_node_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_integer_node_2d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_integer_node_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_integer_node_3d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_integer_node_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_real_cell_1d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_real_cell_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_real_cell_2d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_real_cell_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_real_cell_3d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_real_cell_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_integer_cell_1d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_integer_cell_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_integer_cell_2d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_integer_cell_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_integer_cell_3d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_integer_cell_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functionaldimensionsize(name, dimname, count, ier)
    character(*), intent(in):: name
    character(*), intent(in):: dimname
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_read_grid_functionaldimensionsize_f2c &
      (name, dimname, count, ier)

  end subroutine

  subroutine iricmi_read_grid_functionaldimension_integer(name, dimname, v_arr, ier)
    character(*), intent(in):: name
    character(*), intent(in):: dimname
    integer, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functionaldimension_integer_f2c &
      (name, dimname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functionaldimension_real(name, dimname, v_arr, ier)
    character(*), intent(in):: name
    character(*), intent(in):: dimname
    double precision, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functionaldimension_real_f2c &
      (name, dimname, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functionaltimesize(name, count, ier)
    character(*), intent(in):: name
    integer, intent(out):: count
    integer, intent(out):: ier

    call iricmi_read_grid_functionaltimesize_f2c &
      (name, count, ier)

  end subroutine

  subroutine iricmi_read_grid_functionaltime(name, time_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:), intent(out):: time_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functionaltime_f2c &
      (name, time_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_integer_node_1d(name, dimid, v_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: dimid
    integer, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_integer_node_f2c &
      (name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_integer_node_2d(name, dimid, v_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: dimid
    integer, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_integer_node_f2c &
      (name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_integer_node_3d(name, dimid, v_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: dimid
    integer, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_integer_node_f2c &
      (name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_real_node_1d(name, dimid, v_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: dimid
    double precision, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_real_node_f2c &
      (name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_real_node_2d(name, dimid, v_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: dimid
    double precision, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_real_node_f2c &
      (name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_real_node_3d(name, dimid, v_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: dimid
    double precision, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_real_node_f2c &
      (name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_integer_cell_1d(name, dimid, v_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: dimid
    integer, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_integer_cell_f2c &
      (name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_integer_cell_2d(name, dimid, v_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: dimid
    integer, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_integer_cell_f2c &
      (name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_integer_cell_3d(name, dimid, v_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: dimid
    integer, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_integer_cell_f2c &
      (name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_real_cell_1d(name, dimid, v_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: dimid
    double precision, dimension(:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_real_cell_f2c &
      (name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_real_cell_2d(name, dimid, v_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: dimid
    double precision, dimension(:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_real_cell_f2c &
      (name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_read_grid_functional_real_cell_3d(name, dimid, v_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: dimid
    double precision, dimension(:,:,:), intent(out):: v_arr
    integer, intent(out):: ier

    call iricmi_read_grid_functional_real_cell_f2c &
      (name, dimid, v_arr, ier)

  end subroutine

  subroutine iricmi_write_grid1d_coords(isize, x_arr, ier)
    integer, intent(in):: isize
    double precision, dimension(:), intent(in):: x_arr
    integer, intent(out):: ier

    call iricmi_write_grid1d_coords_f2c &
      (isize, x_arr, ier)

  end subroutine

  subroutine iricmi_write_grid2d_coords_1d(isize, jsize, x_arr, y_arr, ier)
    integer, intent(in):: isize
    integer, intent(in):: jsize
    double precision, dimension(:), intent(in):: x_arr
    double precision, dimension(:), intent(in):: y_arr
    integer, intent(out):: ier

    call iricmi_write_grid2d_coords_f2c &
      (isize, jsize, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_write_grid2d_coords_2d(isize, jsize, x_arr, y_arr, ier)
    integer, intent(in):: isize
    integer, intent(in):: jsize
    double precision, dimension(:,:), intent(in):: x_arr
    double precision, dimension(:,:), intent(in):: y_arr
    integer, intent(out):: ier

    call iricmi_write_grid2d_coords_f2c &
      (isize, jsize, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_write_grid3d_coords_1d(isize, jsize, ksize, x_arr, y_arr, z_arr, ier)
    integer, intent(in):: isize
    integer, intent(in):: jsize
    integer, intent(in):: ksize
    double precision, dimension(:), intent(in):: x_arr
    double precision, dimension(:), intent(in):: y_arr
    double precision, dimension(:), intent(in):: z_arr
    integer, intent(out):: ier

    call iricmi_write_grid3d_coords_f2c &
      (isize, jsize, ksize, x_arr, y_arr, z_arr, ier)

  end subroutine

  subroutine iricmi_write_grid3d_coords_3d(isize, jsize, ksize, x_arr, y_arr, z_arr, ier)
    integer, intent(in):: isize
    integer, intent(in):: jsize
    integer, intent(in):: ksize
    double precision, dimension(:,:,:), intent(in):: x_arr
    double precision, dimension(:,:,:), intent(in):: y_arr
    double precision, dimension(:,:,:), intent(in):: z_arr
    integer, intent(out):: ier

    call iricmi_write_grid3d_coords_f2c &
      (isize, jsize, ksize, x_arr, y_arr, z_arr, ier)

  end subroutine

  subroutine iricmi_write_namedgrid1d_coords(name, isize, x_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: isize
    double precision, dimension(:), intent(in):: x_arr
    integer, intent(out):: ier

    call iricmi_write_namedgrid1d_coords_f2c &
      (name, isize, x_arr, ier)

  end subroutine

  subroutine iricmi_write_namedgrid2d_coords_1d(name, isize, jsize, x_arr, y_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: isize
    integer, intent(in):: jsize
    double precision, dimension(:), intent(in):: x_arr
    double precision, dimension(:), intent(in):: y_arr
    integer, intent(out):: ier

    call iricmi_write_namedgrid2d_coords_f2c &
      (name, isize, jsize, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_write_namedgrid2d_coords_2d(name, isize, jsize, x_arr, y_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: isize
    integer, intent(in):: jsize
    double precision, dimension(:,:), intent(in):: x_arr
    double precision, dimension(:,:), intent(in):: y_arr
    integer, intent(out):: ier

    call iricmi_write_namedgrid2d_coords_f2c &
      (name, isize, jsize, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_write_namedgrid3d_coords_1d(name, isize, jsize, ksize, x_arr, y_arr, z_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: isize
    integer, intent(in):: jsize
    integer, intent(in):: ksize
    double precision, dimension(:), intent(in):: x_arr
    double precision, dimension(:), intent(in):: y_arr
    double precision, dimension(:), intent(in):: z_arr
    integer, intent(out):: ier

    call iricmi_write_namedgrid3d_coords_f2c &
      (name, isize, jsize, ksize, x_arr, y_arr, z_arr, ier)

  end subroutine

  subroutine iricmi_write_namedgrid3d_coords_3d(name, isize, jsize, ksize, x_arr, y_arr, z_arr, ier)
    character(*), intent(in):: name
    integer, intent(in):: isize
    integer, intent(in):: jsize
    integer, intent(in):: ksize
    double precision, dimension(:,:,:), intent(in):: x_arr
    double precision, dimension(:,:,:), intent(in):: y_arr
    double precision, dimension(:,:,:), intent(in):: z_arr
    integer, intent(out):: ier

    call iricmi_write_namedgrid3d_coords_f2c &
      (name, isize, jsize, ksize, x_arr, y_arr, z_arr, ier)

  end subroutine

  subroutine iricmi_read_grid2d_open(grid_handle, ier)
    integer, intent(out):: grid_handle
    integer, intent(out):: ier

    call iricmi_read_grid2d_open_f2c &
      (grid_handle, ier)

  end subroutine

  subroutine iricmi_rin_grid_cell_integer_1d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_cell_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_cell_integer_2d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_cell_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_cell_integer_3d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_cell_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_cell_real_1d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_cell_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_cell_real_2d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_cell_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_cell_real_3d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_cell_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_cell_integer_1d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_cell_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_cell_integer_2d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_cell_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_cell_integer_3d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_cell_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_cell_real_1d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_cell_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_cell_real_2d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_cell_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_cell_real_3d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_cell_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid2d_coords_1d(x_arr, y_arr, ier)
    double precision, dimension(:), intent(in):: x_arr
    double precision, dimension(:), intent(in):: y_arr
    integer, intent(out):: ier

    call iricmi_rin_grid2d_coords_f2c &
      (x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid2d_coords_2d(x_arr, y_arr, ier)
    double precision, dimension(:,:), intent(in):: x_arr
    double precision, dimension(:,:), intent(in):: y_arr
    integer, intent(out):: ier

    call iricmi_rin_grid2d_coords_f2c &
      (x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid3d_coords_1d(x_arr, y_arr, z_arr, ier)
    double precision, dimension(:), intent(in):: x_arr
    double precision, dimension(:), intent(in):: y_arr
    double precision, dimension(:), intent(in):: z_arr
    integer, intent(out):: ier

    call iricmi_rin_grid3d_coords_f2c &
      (x_arr, y_arr, z_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid3d_coords_3d(x_arr, y_arr, z_arr, ier)
    double precision, dimension(:,:,:), intent(in):: x_arr
    double precision, dimension(:,:,:), intent(in):: y_arr
    double precision, dimension(:,:,:), intent(in):: z_arr
    integer, intent(out):: ier

    call iricmi_rin_grid3d_coords_f2c &
      (x_arr, y_arr, z_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid2d_coords_1d(x_arr, y_arr, ier)
    double precision, dimension(:), intent(in):: x_arr
    double precision, dimension(:), intent(in):: y_arr
    integer, intent(out):: ier

    call iricmi_rout_grid2d_coords_f2c &
      (x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid2d_coords_2d(x_arr, y_arr, ier)
    double precision, dimension(:,:), intent(in):: x_arr
    double precision, dimension(:,:), intent(in):: y_arr
    integer, intent(out):: ier

    call iricmi_rout_grid2d_coords_f2c &
      (x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid3d_coords_1d(x_arr, y_arr, z_arr, ier)
    double precision, dimension(:), intent(in):: x_arr
    double precision, dimension(:), intent(in):: y_arr
    double precision, dimension(:), intent(in):: z_arr
    integer, intent(out):: ier

    call iricmi_rout_grid3d_coords_f2c &
      (x_arr, y_arr, z_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid3d_coords_3d(x_arr, y_arr, z_arr, ier)
    double precision, dimension(:,:,:), intent(in):: x_arr
    double precision, dimension(:,:,:), intent(in):: y_arr
    double precision, dimension(:,:,:), intent(in):: z_arr
    integer, intent(out):: ier

    call iricmi_rout_grid3d_coords_f2c &
      (x_arr, y_arr, z_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_iface_integer_1d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_iface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_iface_integer_2d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_iface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_iface_integer_3d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_iface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_iface_real_1d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_iface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_iface_real_2d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_iface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_iface_real_3d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_iface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_iface_integer_1d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_iface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_iface_integer_2d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_iface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_iface_integer_3d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_iface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_iface_real_1d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_iface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_iface_real_2d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_iface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_iface_real_3d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_iface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_jface_integer_1d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_jface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_jface_integer_2d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_jface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_jface_integer_3d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_jface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_jface_real_1d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_jface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_jface_real_2d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_jface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_jface_real_3d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_jface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_jface_integer_1d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_jface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_jface_integer_2d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_jface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_jface_integer_3d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_jface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_jface_real_1d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_jface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_jface_real_2d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_jface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_jface_real_3d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_jface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_kface_integer_1d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_kface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_kface_integer_2d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_kface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_kface_integer_3d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_kface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_kface_real_1d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_kface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_kface_real_2d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_kface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_kface_real_3d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_kface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_kface_integer_1d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_kface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_kface_integer_2d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_kface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_kface_integer_3d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_kface_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_kface_real_1d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_kface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_kface_real_2d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_kface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_kface_real_3d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_kface_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_node_integer_1d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_node_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_node_integer_2d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_node_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_node_integer_3d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_node_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_node_real_1d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_node_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_node_real_2d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_node_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_node_real_3d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_node_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_node_integer_1d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_node_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_node_integer_2d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_node_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_node_integer_3d(name, v_arr, ier)
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_node_integer_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_node_real_1d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_node_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_node_real_2d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_node_real_f2c &
      (name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_node_real_3d(name, v_arr, ier)
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_node_real_f2c &
      (name, v_arr, ier)

  end subroutine



  ! from iricmi_sol_cell.h

  subroutine iricmi_rin_grid_cell_integer_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_cell_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_cell_integer_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_cell_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_cell_integer_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_cell_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_cell_real_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_cell_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_cell_real_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_cell_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_cell_real_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_cell_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_cell_integer_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_cell_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_cell_integer_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_cell_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_cell_integer_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_cell_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_cell_real_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_cell_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_cell_real_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_cell_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_cell_real_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_cell_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine



  ! from iricmi_sol_gridcoord.h

  subroutine iricmi_rin_grid2d_coords_withgridid_1d(gid, x_arr, y_arr, ier)
    integer, intent(in):: gid
    double precision, dimension(:), intent(in):: x_arr
    double precision, dimension(:), intent(in):: y_arr
    integer, intent(out):: ier

    call iricmi_rin_grid2d_coords_withgridid_f2c &
      (gid, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid2d_coords_withgridid_2d(gid, x_arr, y_arr, ier)
    integer, intent(in):: gid
    double precision, dimension(:,:), intent(in):: x_arr
    double precision, dimension(:,:), intent(in):: y_arr
    integer, intent(out):: ier

    call iricmi_rin_grid2d_coords_withgridid_f2c &
      (gid, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid3d_coords_withgridid_1d(gid, x_arr, y_arr, z_arr, ier)
    integer, intent(in):: gid
    double precision, dimension(:), intent(in):: x_arr
    double precision, dimension(:), intent(in):: y_arr
    double precision, dimension(:), intent(in):: z_arr
    integer, intent(out):: ier

    call iricmi_rin_grid3d_coords_withgridid_f2c &
      (gid, x_arr, y_arr, z_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid3d_coords_withgridid_3d(gid, x_arr, y_arr, z_arr, ier)
    integer, intent(in):: gid
    double precision, dimension(:,:,:), intent(in):: x_arr
    double precision, dimension(:,:,:), intent(in):: y_arr
    double precision, dimension(:,:,:), intent(in):: z_arr
    integer, intent(out):: ier

    call iricmi_rin_grid3d_coords_withgridid_f2c &
      (gid, x_arr, y_arr, z_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid2d_coords_withgridid_1d(gid, x_arr, y_arr, ier)
    integer, intent(in):: gid
    double precision, dimension(:), intent(in):: x_arr
    double precision, dimension(:), intent(in):: y_arr
    integer, intent(out):: ier

    call iricmi_rout_grid2d_coords_withgridid_f2c &
      (gid, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid2d_coords_withgridid_2d(gid, x_arr, y_arr, ier)
    integer, intent(in):: gid
    double precision, dimension(:,:), intent(in):: x_arr
    double precision, dimension(:,:), intent(in):: y_arr
    integer, intent(out):: ier

    call iricmi_rout_grid2d_coords_withgridid_f2c &
      (gid, x_arr, y_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid3d_coords_withgridid_1d(gid, x_arr, y_arr, z_arr, ier)
    integer, intent(in):: gid
    double precision, dimension(:), intent(in):: x_arr
    double precision, dimension(:), intent(in):: y_arr
    double precision, dimension(:), intent(in):: z_arr
    integer, intent(out):: ier

    call iricmi_rout_grid3d_coords_withgridid_f2c &
      (gid, x_arr, y_arr, z_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid3d_coords_withgridid_3d(gid, x_arr, y_arr, z_arr, ier)
    integer, intent(in):: gid
    double precision, dimension(:,:,:), intent(in):: x_arr
    double precision, dimension(:,:,:), intent(in):: y_arr
    double precision, dimension(:,:,:), intent(in):: z_arr
    integer, intent(out):: ier

    call iricmi_rout_grid3d_coords_withgridid_f2c &
      (gid, x_arr, y_arr, z_arr, ier)

  end subroutine



  ! from iricmi_sol_iface.h

  subroutine iricmi_rin_grid_iface_integer_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_iface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_iface_integer_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_iface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_iface_integer_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_iface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_iface_real_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_iface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_iface_real_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_iface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_iface_real_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_iface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_iface_integer_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_iface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_iface_integer_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_iface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_iface_integer_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_iface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_iface_real_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_iface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_iface_real_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_iface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_iface_real_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_iface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine



  ! from iricmi_sol_jface.h

  subroutine iricmi_rin_grid_jface_integer_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_jface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_jface_integer_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_jface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_jface_integer_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_jface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_jface_real_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_jface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_jface_real_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_jface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_jface_real_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_jface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_jface_integer_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_jface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_jface_integer_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_jface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_jface_integer_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_jface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_jface_real_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_jface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_jface_real_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_jface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_jface_real_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_jface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine



  ! from iricmi_sol_kface.h

  subroutine iricmi_rin_grid_kface_integer_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_kface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_kface_integer_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_kface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_kface_integer_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_kface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_kface_real_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_kface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_kface_real_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_kface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_kface_real_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_kface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_kface_integer_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_kface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_kface_integer_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_kface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_kface_integer_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_kface_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_kface_real_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_kface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_kface_real_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_kface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_kface_real_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_kface_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine



  ! from iricmi_sol_node.h

  subroutine iricmi_rin_grid_node_integer_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_node_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_node_integer_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_node_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_node_integer_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_node_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_node_real_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_node_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_node_real_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_node_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rin_grid_node_real_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rin_grid_node_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_node_integer_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_node_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_node_integer_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_node_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_node_integer_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    integer, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_node_integer_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_node_real_withgridid_1d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_node_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_node_real_withgridid_2d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_node_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine

  subroutine iricmi_rout_grid_node_real_withgridid_3d(gid, name, v_arr, ier)
    integer, intent(in):: gid
    character(*), intent(in):: name
    double precision, dimension(:,:,:), intent(in):: v_arr
    integer, intent(out):: ier

    call iricmi_rout_grid_node_real_withgridid_f2c &
      (gid, name, v_arr, ier)

  end subroutine



  ! from iricmi_time.h

  subroutine iricmi_rin_time(t, ier)
    double precision, intent(in):: t
    integer, intent(out):: ier

    call iricmi_rin_time_f2c &
      (t, ier)

  end subroutine

  subroutine iricmi_rin_dump_interval(interval, ier)
    double precision, intent(in):: interval
    integer, intent(out):: ier

    call iricmi_rin_dump_interval_f2c &
      (interval, ier)

  end subroutine

  subroutine iricmi_rout_time(t, ier)
    double precision, intent(in):: t
    integer, intent(out):: ier

    call iricmi_rout_time_f2c &
      (t, ier)

  end subroutine

  subroutine iricmi_rout_exchange_interval(interval, ier)
    double precision, intent(in):: interval
    integer, intent(out):: ier

    call iricmi_rout_exchange_interval_f2c &
      (interval, ier)

  end subroutine

  subroutine iricmi_rout_dump_interval(interval, ier)
    double precision, intent(in):: interval
    integer, intent(out):: ier

    call iricmi_rout_dump_interval_f2c &
      (interval, ier)

  end subroutine


end module
