#pragma once

namespace ran_forest 
{
  enum SplittingOrder {BFS,DFS}; // Breadth First Splitting, Depth First Splitting
  enum ElectionStatus { SUCCESS, NODE_SIZE_LIMIT_REACHED, CONVERGED, MAX_DEPTH_REACHED, NULL_HYPOTHESIS_SET };
};

