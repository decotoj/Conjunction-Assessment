# Conjunction-Assessment
All vs All Conjunction Assessment Demo

Uses python's native ability to very quickly sort large lists (in this case lists of positon vector components) to enable fast satellite conjunction assessment with large number of objects.

Example real world TLE catalog included with 20,0451 unique objects.  On the author's desktop class machine a 24 hour all vs all CA run with this catalog and 1 second time steps is currently taking 45 minutes to complete.  

This process will do a leak proof search for conjunctions below a threshold distance (set to 10km by default). As that threshold is increased the process approaches O(n^2).  With the real world current space catalog and a threshold of 10km the process is close to O(n).
