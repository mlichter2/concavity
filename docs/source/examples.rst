.. _examples:

********
EXAMPLES
********

Create a concave hull polygon from a set of coordinates based on k-nearest neighbors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

create a concave hull polygon based on 20 nearest neighbors
-----------------------------------------------------------

create a synthetic set of points::

    import numpy as np
    from shapely.geometry import Point
    import geopandas as gpd
    import matplotlib.pyplot as plt
    import concavity

    np.random.seed(seed=456)
    x = np.random.normal(100000,500, 1000)
    y = np.random.normal(100000, 500, 1000)
    coords = np.array(list(zip(x, y)))
    gpd.GeoSeries([Point(c) for c in coords]).plot()

.. image:: /../figs/points.png
   :align: center
   :alt:

create a concave hull polygon based on 20 nearest neighbors::

   ch = concavity.concave_hull(coords, 20)
   ch




.. image:: /../figs/poly.png
   :align: center
   :alt:


plot concave hull polygon based on 50 nearest neighbors and compare with a convex hull polygon
----------------------------------------------------------------------------------------------
 ::

    ch = concavity.concave_hull(coords, 50)
    concavity.plot_concave_hull(coords, ch)



.. image:: /../figs/knn50.png
   :align: center
   :alt:


create a more elaborate concave hull polygon based on 30 nearest neighbors and compare with a convex hull polygon using plot_concave_hull
-----------------------------------------------------------------------------------------------------------------------------------------
 ::

    ch = concavity.concave_hull(coords, 30)
    concavity.plot_concave_hull(coords, ch)



.. image:: /../figs/knn30.png
   :align: center
   :alt:


create an even more elaborate boundary based on 15 nearest neighbors
--------------------------------------------------------------------
 ::

    ch = concavity.concave_hull(coords, 15)
    concavity.plot_concave_hull(coords, ch)



.. image:: /../figs/knn15.png
   :align: center
   :alt:


Find concave and convex vertices on a polygon boundary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The ``find_concave_vertices`` and ``find_concave_vertices`` functions take a polygon, an angle_threshold as the minimum
angle above which a vertex is in a concave/convex location and a filter type that determines if the function outputs all
 point above the angle threshold or will attempt to locate the peak concave/convex vertices. the output type can be
  either ageopandas GeoDataFrame or a list of vertices and the angle of the two edges they connect

Output all the concave/convex vertices on the polygon boundary
---------------------------------------------------------------
 ::

    concave_df = concavity.find_concave_vertices(ch,0, filter_type ='all')
    convex_df = concavity.find_convex_vertices(ch,0, filter_type ='all')
    concavity.plot_vertices(ch, concave_df, convex_df)



.. image:: /../figs/all_poly_not_smoothed.png
   :align: center
   :alt:


To make a less obvious example - let's smooth our polygon and try again
-----------------------------------------------------------------------
 ::

    from concavity.utils import gaussian_smooth_geom
    geom = gaussian_smooth_geom( ch)
    concave_df = concavity.find_concave_vertices(geom,0, filter_type ='all')
    convex_df = concavity.find_convex_vertices(geom,0, filter_type ='all')
    concavity.plot_vertices(geom, concave_df, convex_df)



.. image:: /../figs/smoothed_all.png
   :align: center
   :alt:


Detect only the peak concave and convex vertices
-------------------------------------------------
 ::

    concave_df = concavity.find_concave_vertices(geom,0, filter_type ='peak')
    convex_df = concavity.find_convex_vertices(geom,0, filter_type ='peak')
    concavity.plot_vertices(geom, concave_df, convex_df)



.. image:: /../figs/peak_thresh0.png
   :align: center
   :alt:


Smooth the angles by using the convovle boolean argument and refine even further the vertices that will be marked as peaks
---------------------------------------------------------------------------------------------------------------------------
 ::

    concave_df = concavity.find_concave_vertices(geom,0, filter_type ='peak', convolve = True)
    convex_df = concavity.find_convex_vertices(geom,0, filter_type ='peak', convolve = True)
    concavity.plot_vertices(geom, concave_df, convex_df)



.. image:: /../figs/smooth_convolve.png
   :align: center
   :alt:


Use the angle threshold argument to limit the angles above which a vertex is considered convex/concave
------------------------------------------------------------------------------------------------------
 ::

    concave_df = concavity.find_concave_vertices(geom,angle_threshold=10, filter_type ='peak')
    convex_df = concavity.find_convex_vertices(geom,angle_threshold=10, filter_type ='peak')
    concavity.plot_vertices(geom, concave_df, convex_df)



.. image:: /../figs/smoothed_thresh10.png
   :align: center
   :alt:





1. Moreira, Adriano & Santos, Maribel. (2007). Concave hull: A k-nearest neighbours approach for the computation of the region occupied by a set of points.. 61-68. 

.. toctree::
   :maxdepth: 2
