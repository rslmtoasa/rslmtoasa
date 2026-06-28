.. _reference/report_module:

=======================================
Report Module
=======================================

Overview
========

The **report module** (``report.f90``) provides a hierarchical data structure for organizing,
collecting, and formatting output information from RS-LMTO-ASA calculations. It implements
a tree-based reporting system that can aggregate statistics, timing information, and
calculation results in a structured, human-readable format.

**Key features:**

- Hierarchical organization of data
- Statistical aggregation (sum, mean, min, max)
- Nested report structures for complex data
- Flexible output formatting
- Performance profiling support
- Integration with timer module

Module Structure
================

Type: report
------------

The ``report`` type implements a tree structure where each node can contain:

- A label (name/description)
- Data values (real or integer arrays)
- Child reports (sub-sections)
- Statistical metadata

**Key members:**

.. code-block:: fortran

   type report
      character(len=sl) :: label              ! Report section name
      real(rp), dimension(:), allocatable :: values  ! Data array
      integer, dimension(:), allocatable :: children_values_index
      type(report), dimension(:), allocatable :: children  ! Sub-reports
      logical :: integer_data                 ! Data type flag
   end type report

Main Procedures
===============

Creation and Building
---------------------

**Constructor:**

.. code-block:: fortran

   type(report) function report_constructor(label)
      character(len=*), intent(in) :: label

Creates a new report with the given label.

**Adding data:**

.. code-block:: fortran

   call my_report%add_value(value)           ! Add single value
   call my_report%add_value(label, value)    ! Add labeled value

**Adding sub-sections:**

.. code-block:: fortran

   call parent_report%add_child(child_report)

Creates hierarchical structure.

Data Retrieval
--------------

**Get values:**

.. code-block:: fortran

   values = my_report%get_values()           ! All values
   values = my_report%get_values(label)      ! Values for sub-section

**Get child reports:**

.. code-block:: fortran

   child = my_report%get_child(label)
   child = my_report%get_valid_child(label)  ! Returns valid child or empty

**Query structure:**

.. code-block:: fortran

   n = my_report%get_num_children()
   n = my_report%get_size_values()
   logical_result = my_report%is_leaf()      ! Has no children?
   logical_result = my_report%has_child(label)

Statistical Operations
----------------------

**Aggregation functions:**

.. code-block:: fortran

   sum_val = my_report%r_sum()               ! Sum of all values
   min_val = my_report%r_min()               ! Minimum value
   max_val = my_report%r_max()               ! Maximum value
   mean_val = my_report%mean()               ! Average value
   perc = my_report%sum_perc()               ! Percentage of parent sum

**Call counting:**

.. code-block:: fortran

   n = my_report%ncalls()                    ! Number of calls to this section
   total = my_report%total_ncalls()          ! Including all children

Output and Formatting
---------------------

**Print report:**

.. code-block:: fortran

   call my_report%print_report(unit, indent_level)

Prints formatted report to file unit with proper indentation.

**Get labels:**

.. code-block:: fortran

   labels = my_report%get_labels()           ! All labels in hierarchy
   parent_label = my_report%get_parent()     ! Parent section name

Use Cases
=========

Performance Profiling
---------------------

Track timing for different code sections:

.. code-block:: fortran

   type(report) :: timing_report
   
   timing_report = report('RS-LMTO Timing')
   
   ! Add timing for main sections
   call timing_report%add_value('SCF Iteration', scf_time)
   call timing_report%add_value('Hamiltonian Build', ham_time)
   call timing_report%add_value('Green Function', gf_time)
   
   ! Create sub-report for detailed breakdown
   scf_report = report('SCF Details')
   call scf_report%add_value('Recursion', recursion_time)
   call scf_report%add_value('DOS Integration', dos_time)
   call timing_report%add_child(scf_report)
   
   ! Print hierarchical timing report
   call timing_report%print_report(6, 0)  ! To stdout

Output:

.. code-block:: text

   RS-LMTO Timing:
     SCF Iteration: 45.2 s (65%)
       SCF Details:
         Recursion: 30.1 s (67%)
         DOS Integration: 15.1 s (33%)
     Hamiltonian Build: 12.3 s (18%)
     Green Function: 11.9 s (17%)

Convergence Tracking
--------------------

Monitor SCF convergence:

.. code-block:: fortran

   type(report) :: scf_report
   
   scf_report = report('SCF Convergence')
   
   do iteration = 1, max_iterations
      ! Do SCF iteration
      call scf_report%add_value('Iteration', real(iteration))
      call scf_report%add_value('ΔQ', charge_diff)
      call scf_report%add_value('Energy', total_energy)
   end do
   
   ! Statistics
   mean_delta_q = scf_report%get_child('ΔQ')%mean()
   final_energy = scf_report%get_child('Energy')%get_values()
   final_energy = final_energy(size(final_energy))

Property Summaries
------------------

Organize calculation results:

.. code-block:: fortran

   type(report) :: results
   
   results = report('Calculation Results')
   
   ! Magnetic properties
   mag_report = report('Magnetic Properties')
   call mag_report%add_value('Total moment', total_moment)
   call mag_report%add_value('Spin moment', spin_moment)
   call mag_report%add_value('Orbital moment', orbital_moment)
   call results%add_child(mag_report)
   
   ! Electronic properties
   elec_report = report('Electronic Properties')
   call elec_report%add_value('Fermi energy', fermi_energy)
   call elec_report%add_value('DOS at EF', dos_at_ef)
   call results%add_child(elec_report)

Hierarchical Data Structure
============================

**Tree representation:**

.. code-block:: text

   Root Report
   ├── Child 1
   │   ├── Grandchild 1.1
   │   └── Grandchild 1.2
   ├── Child 2
   └── Child 3
       ├── Grandchild 3.1
       ├── Grandchild 3.2
       └── Grandchild 3.3

**Navigation:**

- Traverse from root to leaves
- Query by label at any level
- Aggregate statistics up the tree
- Parent-child relationships preserved

Implementation Details
======================

**Memory management:**

- Allocatable arrays for flexible sizing
- Dynamic child array growth
- Efficient value storage

**Data types:**

- Real (double precision) values
- Integer values (with flag)
- Mixed types via generic interfaces

**Generic interfaces:**

.. code-block:: fortran

   generic :: add_value => report_add_value, &
                           report_add_label_value, &
                           report_add_value_r4, &
                           report_add_value_i, &
                           ...

Integration with RS-LMTO
=========================

**Timer integration:**

The report module works closely with ``timer_mod`` for performance profiling:

.. code-block:: fortran

   use timer_mod, only: g_timer
   use report_mod
   
   call g_timer%start('SCF Loop')
   ! ... SCF calculation ...
   call g_timer%stop('SCF Loop')
   
   ! Export to report
   timing_report = g_timer%get_report()

**Logger integration:**

Works with ``logger_mod`` for structured output:

.. code-block:: fortran

   use logger_mod, only: g_logger
   
   call results_report%print_report(g_logger%unit, indent=2)

Example: Complete Reporting Workflow
=====================================

.. code-block:: fortran

   program report_example
      use report_mod
      use precision_mod
      implicit none
      
      type(report) :: main_report, scf_report, dos_report
      integer :: iter
      real(rp) :: charge_diff, energy
      
      ! Create main report
      main_report = report('RS-LMTO Calculation')
      
      ! SCF section
      scf_report = report('Self-Consistent Field')
      do iter = 1, 10
         ! Simulate SCF
         charge_diff = 0.1_rp / iter
         energy = -100.0_rp + 1.0_rp / iter
         
         call scf_report%add_value('ΔQ', charge_diff)
         call scf_report%add_value('Energy', energy)
      end do
      call main_report%add_child(scf_report)
      
      ! DOS section
      dos_report = report('Density of States')
      call dos_report%add_value('DOS at EF', 2.3_rp)
      call dos_report%add_value('Integrated DOS', 24.0_rp)
      call main_report%add_child(dos_report)
      
      ! Print full report
      call main_report%print_report(6, 0)
      
      ! Extract specific information
      print *, 'Final ΔQ:', scf_report%get_child('ΔQ')%mean()
      print *, 'Converged energy:', energy
      print *, 'Total items:', main_report%total_ncalls()
      
   end program report_example

Output Format
=============

**Standard format:**

.. code-block:: text

   RS-LMTO Calculation:
     Self-Consistent Field:
       ΔQ: 0.0100 (mean over 10 iterations)
       Energy: -99.10 eV (final value)
     Density of States:
       DOS at EF: 2.30 states/eV
       Integrated DOS: 24.00 electrons

**With percentages:**

.. code-block:: text

   Timing Report:
     Total time: 125.5 s (100%)
       SCF: 89.3 s (71%)
         Recursion: 45.2 s (51% of SCF, 36% of total)
         Green's function: 30.1 s (34% of SCF, 24% of total)
         Mixing: 14.0 s (16% of SCF, 11% of total)
       Post-processing: 36.2 s (29%)

Advanced Features
=================

**Filtering:**

Extract specific data based on labels:

.. code-block:: fortran

   scf_times = main_report%get_child('Timing')%get_child('SCF')%get_values()

**Comparison:**

Compare reports from different calculations:

.. code-block:: fortran

   calc1_energy = report1%get_child('Energy')%get_values()
   calc2_energy = report2%get_child('Energy')%get_values()
   energy_diff = calc2_energy - calc1_energy

**Serialization:**

Reports can be written to files for later analysis:

.. code-block:: fortran

   open(unit=10, file='report.dat')
   call my_report%print_report(10, 0)
   close(10)

Best Practices
==============

1. **Consistent naming:** Use clear, descriptive labels
2. **Hierarchical organization:** Group related data logically
3. **Units in labels:** Include units in labels (e.g., 'Time (s)', 'Energy (eV)')
4. **Regular updates:** Add values incrementally during calculation
5. **Summary at end:** Print full report at calculation completion
6. **Combine with timing:** Use report module with timer for comprehensive profiling

Provenance
==========

**Implementation:**

- ``source/report.f90`` - Main report module (1100+ lines)
- ``source/timer_mod.f90`` - Timing integration
- ``source/logger_mod.f90`` - Output integration
- ``source/string_mod.f90`` - String utilities for formatting

**Dependencies:**

- ``precision_mod`` - Floating-point precision
- ``array_mod`` - Dynamic array operations
- ``string_mod`` - String manipulation

See Also
========

- :ref:`reference/module_overview` - Other utility modules
- ``source/logger.f90`` - Logging infrastructure
- ``source/timer.f90`` - Performance timing
