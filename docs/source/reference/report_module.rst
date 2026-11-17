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

   call my_report%get_values(values)         ! All values (subroutine)
   call my_report%get_values(label, values)  ! Values for sub-section

**Get child reports:**

.. code-block:: fortran

   type(report), pointer :: child
   found = my_report%get_child(label, child)        ! logical result; child is set
   found = my_report%get_valid_child(label, child)  ! logical result; child is set

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

   call my_report%print_report(title)
   call my_report%print_report(title, fields)

Prints a table to standard output. ``fields`` is a comma-separated list;
the default is ``'MEAN,MAX,MIN,TOTAL,TOTAL (%),NCALLS,TOTAL NCALLS'``.

**Get labels:**

.. code-block:: fortran

   call my_report%get_labels(labels)                ! All labels in hierarchy (subroutine)
   found = my_report%get_parent(label, parent)      ! logical result; parent is set

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

``timer_mod`` (``g_timer``) has ``start``, ``stop`` and ``print_report``:

.. code-block:: fortran

   use timer_mod, only: g_timer
   
   call g_timer%start('SCF Loop')
   ! ... SCF calculation ...
   call g_timer%stop('SCF Loop')
   
   call g_timer%print_report()

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

- ``source/report.f90`` - Main report module
- ``source/timer.f90`` - Timing integration
- ``source/logger.f90`` - Output integration
- ``source/string.f90`` - String utilities for formatting

**Dependencies:**

- ``precision_mod`` - Floating-point precision
- ``array_mod`` - Dynamic array operations
- ``string_mod`` - String manipulation

See Also
========

- :ref:`reference/module_overview` - Other utility modules
- ``source/logger.f90`` - Logging infrastructure
- ``source/timer.f90`` - Performance timing
