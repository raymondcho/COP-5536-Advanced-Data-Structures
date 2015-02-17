import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Random;
import java.util.StringTokenizer;

public class dijkstra
{
	/*
	 *  Application settings. 
	 */
	boolean allowimprove = true;			// Toggle this to true if wish to activate improved random edge generator (see documentation).
	boolean outputadjacentlist = false;		// Toggle this to print out adjacency list after graph is generated (for both random and user input modes).
	boolean outputfullresultpath = false;	// Toggle this to print out full path from source vertex to each other vertex.
	boolean suppressresultoutput = false;	// Toggle this to suppress output of dijkstra's algorithm, regardless of outputfullresultpath setting.
	boolean outputresultonly = true;		// Toggle this to output the result of each row only. All other output is suppressed.
	
	class userInputEdges		// Helper data structure for temporary storage of edges from user input mode
	{
		int v1;
		int v2;
		int c;
		public userInputEdges(int v1, int v2, int c)
		{
			this.v1 = v1;
			this.v2 = v2;
			this.c = c;
		}
	}
	
	class Vertex				// Vertex structure holds index and distance values.
	{
		int index;
		int distance;
		public Vertex(int index1, int distance1)
		{
			index = index1;
			distance = distance1;
		}
	}
	
	class Graph					// Graph structure for holding information of the graph.
	{
		class row				// Row structure, for holding each vertex's edge connections to other vertices.
		{

			LinkedList<Vertex> vertices;			// Stores adjacent list.
			HashSet<Integer> connectedgraphs;		// Used for quick check if vertex is already connected to second vertex.
			
			public row()
			{
				vertices = new LinkedList<Vertex>();
				connectedgraphs = new HashSet<Integer>();
			}
			
			public boolean checkedge(int j)
			{
				if (connectedgraphs.contains(j))
				{
					return true;
				}
				else
				{
					return false;
				}
			}
			
			public void addedge(int j, int cost)
			{
				if (!connectedgraphs.contains(Integer.valueOf(j)))
				{
					Vertex a = new Vertex(j, cost);
					vertices.add(a);
					connectedgraphs.add(j);
				}
			}
		}
		
		int n;
		double d;
		int x;
		int edges;
		row[] graphs;
		HashSet<Integer> visited2 = new HashSet<Integer>();
		int mode;
		boolean verified;		// This ensures we don't run either simple or f-heap scheme unless we know graph is fully connected.
		
		public Graph(int n1, double d1, int x1, int m1)
		{
			n = n1;
			d = d1;
			x = x1;
			verified = false;
			mode = m1;
			double tedges = (0.01)*(((n-1)*n)/2);
			edges = (int) (tedges*d);
			graphs = new row[n];
			for (int z = 0; z < n; z++)
			{
				graphs[z] = new row();
			}
			if (!outputresultonly)
			{
				System.out.println("Initializing new graph with "+n+" nodes, edge density of "+d+"% and # of edges = "+edges);
			}
			long st = System.currentTimeMillis();
			if (d < 1)
			{
				if (!outputresultonly)
				{
					System.out.println("Due to extremely low number of edges needed, switching to special random edge generator (see documentation).");
				}
				initializeSparseGraph();
			}
			else if (d == 100 || allowimprove && d > 20)
			{
				if (!outputresultonly)
				{
					System.out.println("Due to high number of edges needed, switching to improved random edge generator (see documentation).");
				}
				improvedInitializeGraph();
			}
			else
			{
				if (!outputresultonly)
				{
					System.out.println("Using standard random edge generator.");
				}
				initializegraph();
			}
			long et = System.currentTimeMillis();
			if (!outputresultonly)
			{
				System.out.println("Finished generating random graph. Graph generation took "+(et-st)+" milliseconds.");
				System.out.println("Checking for fully connected graph ...");
			}
			boolean check = false;
			while (check == false)
			{
				if (initializecheckgraph())
				{
					if (!outputresultonly)
					{
						System.out.println("Graph is fully connected.");
					}
					verified = true;
					check = true;
					if (outputadjacentlist && !outputresultonly)
					{
						System.out.println("Graph adjacency list:");
						int county = 0;
							for (int i = 0; i < n; i++)
							{
								System.out.print("Vertex "+i+":\t");
								ListIterator<Vertex> za = graphs[i].vertices.listIterator();
								while(za.hasNext())
								{
									Vertex nex = za.next();
									System.out.print(nex.index+" ("+nex.distance+")  ");
									county++;
								}
								System.out.println("");
							}
						System.out.println("Total of "+county+" edge pairs");
					}
				}
				else
				{
					if (!outputresultonly)
					{
						System.out.println("Graph is not fully connected!!");
						System.out.println("Re-generating another graph with same input parameters ...");
					}

					graphs = new row[n];
					for (int z = 0; z < n; z++)
					{
						graphs[z] = new row();
					}
					if (!outputresultonly)
					{
						System.out.println("Initializing new graph with "+n+" nodes, edge density of "+d+"% and # of edges = "+edges);
					}
					long st1 = System.currentTimeMillis();
					if (d < 1)
					{
						if (!outputresultonly)
						{
							System.out.println("Due to extremely low number of edges needed, switching to special random edge generator (see documentation).");
						}
						initializeSparseGraph();
					}
					else if (d == 100 || allowimprove && d > 20)
					{
						if (!outputresultonly)
						{
							System.out.println("Due to high number of edges needed, switching to improved random edge generator (see documentation).");
						}
						improvedInitializeGraph();
					}
					else
					{
						if (!outputresultonly)
						{
							System.out.println("Using standard random edge generator.");
						}
						initializegraph();
					}
					long et1 = System.currentTimeMillis();
					if (!outputresultonly)
					{
						System.out.println("Finished generating random graph. Graph generation took "+(et1-st1)+" milliseconds.");
						System.out.println("Checking for fully connected graph ...");
					}
				}
			}
		}
		
		public Graph(int x, int n, int m, userInputEdges[] edges, int m1)
		{
			d = 0;
			this.n = n;
			this.x = x;
			this.edges = m;
			mode = m1;
			graphs = new row[n];
			verified = false;
			for (int z = 0; z < n; z++)
			{
				graphs[z] = new row();
			}
			for (int i = 0; i < edges.length; i++)
			{
				graphs[edges[i].v1].addedge(edges[i].v2, edges[i].c);
				graphs[edges[i].v2].addedge(edges[i].v1, edges[i].c);
			}
			if (!outputresultonly)
			{
				System.out.println("Successfully created user-provided graph.");
				System.out.println("Checking for fully connected graph ...");
			}
			if (initializecheckgraph())
			{
				if (!outputresultonly)
				{
					System.out.println("Graph is fully connected.");
				}
				verified = true;
			}
			else
			{

				System.out.println("Graph is not fully connected!! Please input a graph that is fully connected.");
			}
		}
		
		private void initializegraph()
		{
			/*
			 * The simple edge generator pulls two random vertices and a random cost, checks if the pair
			 * exists in the graph, and if not, adds that pair to the graph. It is not recommended to use
			 * this edge generator in cases where edge density > 20%. Performance begins to noticeably 
			 * degrade as "collisions" -- (re-creating already paired vertex edges) start to take over.
			 * 
			 * The simple edge generator is excellent, however, in cases of low density. In which case
			 * the probability of a collision is minimized since only a few edges (relatively) are created.
			 * However, in cases of very low density, another edge generator should be used, as this one
			 * is not well suited in cases of generating a fully connected sparse graph. 
			 * See initializeSparseGraph() below.
			 */
			int counter = 0;											// Counts number of edges created.
			int initialedgesrange = edges-counter;
			if (!outputresultonly)
			{
				System.out.print("Status: ");
				System.out.println("------------------------|25%---------------------|50%---------------------|75%---------------------|100%");
				System.out.print("        ");
			}

			int increment = initialedgesrange/100;
			int remainder = initialedgesrange % 100;
			int runningincrement = increment+remainder;
			while (counter < edges)
			{
				int i = randomi(n);									// Random value for v1.
				int j = randomi(n);									// Random value for v2.
				int cost = randomi(1000)+1;							// Random value for cost/distance.
				
				if (!graphs[i].checkedge(j) && !graphs[j].checkedge(i))	// Only add edge if it doesn't already exist.
				{
					graphs[i].addedge(j, cost);
					graphs[j].addedge(i, cost);
					if (counter > runningincrement && !outputresultonly)
					{
						System.out.print(".");
						runningincrement = runningincrement+increment;
					}
					counter++;
				}
			}
			if (!outputresultonly)
			{
				System.out.println("Done!");
			}
			
		}
		
		private void improvedInitializeGraph()
		{
			/*
			 * The improved edge generator is especially adapted to handle high density edge cases.
			 * The main drawback for the standard edge generator, is that when edge density gets higher
			 * it takes an exponentially longer time to generate new edges because the probability that
			 * the same edge pair (either between v1 and v2 or v2 and v1) has already been generated 
			 * increases exponentially. 
			 * Thus, the improved edge generator starts from a list of *all* possible pair connections
			 * each with a random cost. Then, if density is 100%, then we are done. If density is < 100%
			 * then we calculate the number of edges that need to be removed, then randomly remove the 
			 * edges from the list and then output the list of edges to the graph generator.
			 * This drastically speeds up random graph generation for high density cases, while still
			 * maintaining a "random" edge generating scheme.
			 * However, in cases of low edge density, the simple/standard edge generator should be used,
			 * to avoid having to compute all pairs only to then take out a vast majority of them, to
			 * leave a sparse graph.
			 */
			ArrayList<userInputEdges> ue = new ArrayList<userInputEdges>();
			for (int i = 0; i < n; i++)
			{
				for (int j = i+1; j < n; j++)
				{
					if (i != j)
					{
						int rc = randomi(1000)+1;							// Random value for cost/distance.
						ue.add(new userInputEdges(i, j, rc));
					}
				}
			}
			int redges = (int) ((((n-1)*n)/2));
			int difference = redges - edges;
			int counter = 0;
			int initialedgesrange = difference;
			if (!outputresultonly)
			{
				System.out.print("Status: ");
				System.out.println("------------------------|25%---------------------|50%---------------------|75%---------------------|100%");
				System.out.print("        .");
			}

			if (difference > 0)
			{
				int increment = initialedgesrange/100;
				int remainder = initialedgesrange % 100;
				int runningincrement = increment+remainder;
				while (counter < difference && !ue.isEmpty())
				{
					int range = ue.size();
					int roll = randomi(range);
					ue.set(roll, null);
					counter++;
					if (counter > runningincrement && !outputresultonly)
					{
						System.out.print(".");
						runningincrement = runningincrement+increment;
					}
				}
				if (!outputresultonly)
				{
					System.out.println("Done! Now outputting edges to graph generator.");
				}
				int increm = ue.size()/100;
				int incremr = ue.size() % 100;
				int runningincrem = increm + incremr;
				if (!outputresultonly)
				{
					System.out.print("Status: ");
					System.out.println("------------------------|25%---------------------|50%---------------------|75%---------------------|100%");
					System.out.print("        .");
				}

				for (int i = 0; i < ue.size(); i++)
				{
					userInputEdges aq = ue.get(i);
					if (aq != null)
					{
						int tv1 = aq.v1;
						int tv2 = aq.v2;
						int tc = aq.c;
						graphs[tv1].addedge(tv2, tc);
						graphs[tv2].addedge(tv1, tc);
					}
					if (i > runningincrem && !outputresultonly)
					{
						System.out.print(".");
						runningincrem = runningincrem+increm;
					}
				}
				if (!outputresultonly)
				{
					System.out.println("Done!");
				}
			}
			else
			{
				if (!outputresultonly)
				{
					System.out.print("............................................................");
					System.out.println("......................................Done! Now outputting edges to graph generator.");
					System.out.print("Status: ");
					System.out.println("------------------------|25%---------------------|50%---------------------|75%---------------------|100%");
					System.out.print("        .");
				}

				int increm = ue.size()/100;
				int incremr = ue.size() % 100;
				int runningincrem = increm + incremr;
				for (int i = 0; i < ue.size(); i++)
				{
					userInputEdges aq = ue.get(i);
					int tv1 = aq.v1;
					int tv2 = aq.v2;
					int tc = aq.c;
					graphs[tv1].addedge(tv2, tc);
					graphs[tv2].addedge(tv1, tc);
					if (i > runningincrem && !outputresultonly)
					{
						System.out.print(".");
						runningincrem = runningincrem+increm;
					}
				}
				if (!outputresultonly)
				{
					System.out.println("Done!");
				}
			}
		}

		private void initializeSparseGraph()
		{
			/*
			 * This edge generator is especially suited for generating graphs with < 1% of edge density.
			 * With edge density < 1%, it becomes very unlikely the standard edge generator will generate
			 * a fully connected graph. 
			 * 
			 * At the initial pass, this generator pairs up each vertex with a random partner vertex,
			 * and assigns a random cost/distance to the pair. Then, for the first vertex, that vertex 
			 * is added to a second pile. That second pile is randomized and two pairs are then pulled
			 * out and matched with a random cost. The first of that pair is put pack into the pile again.
			 * This is done recursively until the density ratio is met.
			 * 
			 * Of course, using this method will eventually exhaust the pile since each successive pile
			 * is only half the size of the previous pile. More specifically, we can get about slightly
			 * less than n number of edges out of this method (which will usually still fall less than
			 * the number of edges to meet density requirement). That is why when the pile becomes < 2,
			 * we switch to the default standard random edge generation to fill in the rest.
			 */
			int counter = 0;											// Counts number of edges created.
			int initialedgesrange = edges-counter;
			if (!outputresultonly)
			{
				System.out.print("Status: ");
				System.out.println("------------------------|25%---------------------|50%---------------------|75%---------------------|100%");
				System.out.print("        ");
			}
			
			int increment = initialedgesrange/100;
			int remainder = initialedgesrange % 100;
			int runningincrement = increment+remainder;
			ArrayList<Integer> pile1 = new ArrayList<Integer>();
			for (int i = 0; i < n; i++)
			{
				pile1.add(Integer.valueOf(i));					// Add all node indices into pile.
			}
			Collections.shuffle(pile1);							// Shuffle the order of the indices.
			ArrayList<Integer> pile2 = null;
			while (counter < edges)
			{
				pile2 = new ArrayList<Integer>();
				while (pile1.size() >= 2 && counter < edges)
				{
					int i = pile1.remove(pile1.size()-1);		// Remove the last two indices.
					int j = pile1.remove(pile1.size()-1);		// Pair them up with a random cost/distance.
					int cost = randomi(1000)+1;	
					pile2.add(Integer.valueOf(i));
					graphs[i].addedge(j, cost);					// Add the pair of edges into the graph.
					graphs[j].addedge(i, cost);
					
					if (counter > runningincrement && !outputresultonly)
					{
						System.out.print(".");
						runningincrement = runningincrement+increment;
					}
					counter++;
				}
				if (!pile1.isEmpty() && counter < edges)
				{
					int o1 = pile1.remove(pile1.size()-1);
					pile2.add(Integer.valueOf(o1));
					int o2 = randomi(n);
					while (o2 == o1)			// Make sure the vertices are not the same.
					{
						o2 = randomi(n);
					}
					int c2 = randomi(1000)+1;	// Generate random cost/distance.
					graphs[o1].addedge(o2, c2);	// Add pair of edges to graph.
					graphs[o2].addedge(o1, c2);
					if (counter > runningincrement && !outputresultonly)
					{
						System.out.print(".");
						runningincrement = runningincrement+increment;
					}
					counter++;
				}
				if (pile2.size() > 1 && counter < edges)
				{
					pile1 = new ArrayList<Integer>(pile2.size());
					for (int yz = 0; yz < pile2.size(); yz++)
					{
						Integer as = pile2.get(yz);
						pile1.add(as);
					}
					pile2.clear();
					Collections.shuffle(pile1);
				}
				
				// Once second pile has been reduced to 1, switch to default standard edge generation. 
				int i = randomi(n);									// Random value for v1.
				int j = randomi(n);									// Random value for v2.
				int cost = randomi(1000)+1;							// Random value for cost/distance.
				
				if (!graphs[i].checkedge(j) && !graphs[j].checkedge(i))	// Only add edge if it doesn't already exist.
				{
					graphs[i].addedge(j, cost);
					graphs[j].addedge(i, cost);
					if (counter > runningincrement && !outputresultonly)
					{
						System.out.print(".");
						runningincrement = runningincrement+increment;
					}
					counter++;
				}
			}
			if (!outputresultonly)
			{
				System.out.println("Done!");
			}
			
		}
		
		public boolean initializecheckgraph()		// Returns true if all nodes can be reached from the source node, false otherwise.
		{
			visited2.clear();
			for (int i = 0; i < n; i++)				// Populate set with all node numbers.
			{
				visited2.add(Integer.valueOf(i));
			}
			checkgraphconnected(x);					// Start depth-first search from the source node.
			if (visited2.isEmpty())					// If there are still nodes remaining, then those nodes are isolated and cannot
			{										// be reached from the source node. Graph is not fully connected!
				return true;
			}
			return false;
		}
		
		public void checkgraphconnected(int ind)	// Helper method uses Depth-First-Search to traverse graph.
		{
			if (visited2.isEmpty() || !visited2.contains(Integer.valueOf(ind)))		// If either condition is encountered, skip this path.
			{
				return;
			}
			else									// For each unvisited node, remove it from set of unvisited nodes, 
			{										// then check all of its neighbors.
				visited2.remove(Integer.valueOf(ind));
				ListIterator<Vertex> li2 = graphs[ind].vertices.listIterator();
				while (li2.hasNext())
				{
					checkgraphconnected(li2.next().index);
				}
			}
		}
	}
	
	class simpleScheme			// Uses simple data structures for O(n^2) runtime
	{
		class simplePriorityQ		// Simple priority queue using a static array for determining next closest vertex.
		{							// This helper method extracts and returns the vertex in Q with smallest cost in cost[].
			int[] backend;
			int size;
			public simplePriorityQ(int size)
			{
				this.size = size;
				backend = new int[size];
				for (int i = 0; i < size; i++)
				{
					backend[i] = 1;
				}
			}
			public boolean isEmpty()
			{
				return (size == 0);
			}
			public int getNextClosest()		// Each call involves a scan of entire array for the next closest vertex.
			{
				if (size > 0)
				{
					int closest = Integer.MAX_VALUE;
					int closestvertex = -1;
					for (int i = 0; i < backend.length; i++)
					{
						if (backend[i] == 1 && cost[i] < closest)		// Only indexes with value = 1 are considered.
						{
							closest = cost[i];
							closestvertex = i;
						}
					}
					backend[closestvertex] = 0;		// The closest vertex is "removed" by setting its value to 0.
					Qa.remove(Integer.valueOf(closestvertex));
					size--;
					return closestvertex;
				}
				else
				{
					return -1;		// Means there are no more vertices to return.
				}
			}
		}
		
		Graph igraph;				// The input graph
		int startindex;				// Index to the source node
		int[] cost;					// Array contains each node's cost/distance for dijkstra's algorithm.
		int[] previous;				// Each node's previous node for dijkstra's algorithm.
		simplePriorityQ Q;			// Simple priority queue using a static array.
		HashSet<Integer> Qa;		// Helper data structure for quick "contains" checking.
		long runningtime;			
		
		public simpleScheme(Graph graph)	// Constructor
		{
			igraph = graph;
			startindex = graph.x;
			cost = new int[graph.n];
			previous = new int[graph.n];
			Qa = new HashSet<Integer>();
			runningtime = 0;
			
			if (!graph.verified)
			{
				if (!outputresultonly)
				{
					System.out.println("Graph has not been verified to be fully connected, or it has been determined it is not fully connected. Abort.");
				}
			}
			else
			{
				executeDijkstra();			// Starts the Dijkstra algorithm.
			}
		}
		
		public void executeDijkstra()
		{
			long starttime = System.currentTimeMillis();		// Starts the measurement of execution time.
			
			for (int i = 0; i < igraph.n; i++)	// Initialization step. Set all costs to infinity and previous nodes to undefined.
			{
				cost[i] = Integer.MAX_VALUE;	// Highest integer value chosen to mean "infinity".
				previous[i] = -99;				// Arbitrary value < 0 chosen to mean "undefined".
			}
			
			cost[startindex] = 0;				// Set the distance from source to source as 0.
			Q = new simplePriorityQ(igraph.n);	// At start all nodes unoptimized, so all are in Q.
			for (int i = 0; i < igraph.n; i++)	// At start all nodes are also in helper set Qa.
			{
				Qa.add(i);
			}
			
			while (!Q.isEmpty())				// The main loop
			{
				int u = Q.getNextClosest();		// Calls helper method to get the next closest node to work on.
				
				if (cost[u] == Integer.MAX_VALUE)
				{
					break;						// Means node is unreachable, end algorithm.
				}
				
				LinkedList<Vertex> V = igraph.graphs[u].vertices;		// Get all neighbors of node u
				ListIterator<Vertex> l2a = V.listIterator();
				
				while (l2a.hasNext())		// For each neighbor v of u where v has not been removed from Q.
				{
					Vertex v = l2a.next();
					if (Qa.contains(Integer.valueOf(v.index)))	// Checks if node v is still in Q via helper Qa.
					{
						int temp = cost[u] + v.distance;		// Relax (u,v).
						if (temp < cost[v.index])				// Update cost/distance if necessary.
						{
							cost[v.index] = temp;
							previous[v.index] = u;
						}
					}
				}
				
			}
			
			long endtime = System.currentTimeMillis();
			runningtime = endtime - starttime;
			
			// Output results.
			if (outputfullresultpath && !outputresultonly)
			{
				for (int i = 0; i < igraph.n; i++)			
				{
					System.out.print(cost[i]+"\t// cost from node "+igraph.x+" to "+i);
					int backtrack = previous[i];
					while (backtrack != -99)
					{
						System.out.print(" from "+backtrack);
						backtrack = previous[backtrack];
					}
					System.out.println("");
				}
			}
			else if (suppressresultoutput && !outputresultonly)
			{
				
			}
			else if (outputresultonly)
			{
				for (int i = 0; i < igraph.n; i++)
				{
					System.out.println(cost[i]);
				}
			}
			else
			{
				for (int i = 0; i < igraph.n; i++)			
				{
					System.out.println(cost[i]+"\t// cost from node "+igraph.x+" to "+i);
				}
			}
			
		}
	}
	
	class doubleFibChain	// List data structure to hold nodes. Has O(1) amortized add, remove, get.
	{
		ArrayList<fibNode> backend;		// ArrayList as underlying data structure.
		int size;				// Modifications have been made to ensure O(1) amortized add, remove, and get.
		int bindex;				// See below for explanation.
		
		public doubleFibChain()
		{
			size = 0;
			bindex = 0;
			backend = new ArrayList<fibNode>();
		}
		
		public void add(fibNode f)
		{
			backend.add(f);					// Adds the node to the list. This has amortized O(1) cost.
			f.fibChainindex = bindex++;		// Sets the node's chain index to its index in the list.
			size++;
		}
		
		public fibNode removeLast()			// Removes and returns the last element of the list.
		{
			/*
			 * Because this method is only called when we are iterating through the entire list of the root or a child list,
			 * we use this time to finish removing any and all nodes that have been nullified by the remove() method.
			 * Essentially, while removing and returning the last node of the list, we will remove, but won't return, nulled entries.
			 * 
			 * Removing from the end of the ArrayList ensures that no objects are shifted down, which is costly. 
			 * Removing from the end of the ArrayList is O(1) cost.
			 */
			fibNode result = null;			// Once the chain is empty, it will return a null fibNode, which signals chain is empty.
			while (result == null && !backend.isEmpty())
			{
				result = backend.remove(backend.size()-1);
			}
			bindex = backend.size();
			return result;
		}
		
		
		public void remove(fibNode f)
		{	
			/*
			 * 	This method "removes" the node from the list (root or child) by nulling its index on the list.
			 *  This ensures the on-spot "removal" step is O(1).
			 *  We only null out the value instead of removing the object from the ArrayList in order to avoid
			 *  having to shift objects from their index position, which is costly. The object can be nulled out
			 *  in O(1) time by calling the fibNode's index field.
			 *  Then when the list (root or child) is iterated over either in the consolidation step for the root
			 *  list or when the parent node is extracted and the children are put onto the root list, null entries
			 *  are skipped over at that point. It matters little that there is an extra few null values when iterating
			 *  through the list. But it matters much that we do as little work as possible at the point of doing
			 *  the on-spot "removal" of the node from the list.
			 */
			backend.set(f.fibChainindex, null);		// "Removes" the node from the list by nulling its index on the list.
													// This ensures the on-spot "removal" step is O(1).
			size--;									// Decrement size of f-heap.
		}
		
		public void clear()		// Empties out the root list.
		{
			backend.clear();
			bindex = backend.size();
			size = 0;
		}
	}
	
	class fibNode			// Class for the fibNode structure. 
	{
		int fibChainindex;	// This node's index in the list (root list or child list of a parent node).
		fibNode parent;		// Pointer to node's parent.
		doubleFibChain child;	// Child list for this node's children.
		boolean mark;		// Mark for purposes of the Fibonacci heap.
		int key;			// Designates the node's unique index number.
		int value;			// The node's value or cost or distance.
		int degree;			// Node's degree = # of immediate children.
		
		public fibNode(int key, int val)
		{
			parent = null;
			mark = false;
			this.key = key;
			value = val;
			degree = 0;
			child = null;
			fibChainindex = 0;
		}
	}
	
	class fibHeap			// Class for the f-heap structure.
	{
		fibNode min;		// Pointer to minimum fibNode.
		int size;			// Stores current size of heap.
		doubleFibChain roots;	// Root list of heap.
		
		public fibHeap()
		{
			min = null;
			size = 0;
			roots = new doubleFibChain();
		}
	}
	
	class FibonacciHeap				// Overall class for implementing the Fibonacci Heap structure
	{
		public FibonacciHeap()
		{
		}
		
		public fibHeap Make_Fib_Heap()
		{
			fibHeap a = new fibHeap();						// Creates and returns a new empty f-heap.
			return a;
		}

		public void Fib_Heap_Insert(fibHeap H, fibNode x)
		{
			H.roots.add(x);									// Inserts the new fibNode into the root list.
			if (H.min == null || x.value < H.min.value)		// If the new node's value is < current min, replace min.
			{
				H.min = x;
			}
			H.size++;										// Increase size of f-heap by 1.
		}
		
		public fibNode Minimum(fibHeap H)
		{
			return H.min;									// Returns the current minimum fibNode.
		}
		
		public fibNode Extract_Min(fibHeap H)				// Pulls and returns the minimum fibNode and consolidates root list.
		{
			fibNode z = H.min;
			if (z != null)
			{
				// Remove minimum node from the root list.
				H.roots.remove(z);

				// Check if old MinRoot node has children node(s). If so, put them all into root list.
				if (z.child != null)
				{
					fibNode c = z.child.removeLast();
					while (c != null)
					{
						c.parent = null;					// Strip each child node of its parent pointer.
						H.roots.add(c);
						c = z.child.removeLast();
					}
				}

				Consolidate(H);								// Consolidate the root list.
				
				H.size--;									// Decrement the f-heap size.
			}
			return z;
		}
		
		public void Consolidate(fibHeap H)		// Consolidation step to merge trees of same degree.
		{
			fibNode[] A = new fibNode[30];		// Temporary array for consolidation step.
			for (int z = 0; z < A.length; z++)
			{
				A[z] = null;
			}

			fibNode x = H.roots.removeLast();	// Start deconstructing root list and putting each tree into temporary array.
			fibNode y = null;
	
			while (x != null)
			{
				int d = x.degree;				// Obtain the node's degree.
				while (A[d] != null)			// This means we already have a tree of same degree. So merge the two trees.
				{
					y = A[d];
					if (x.value > y.value)		// If y has a lesser value than x, swap the two before calling linker method.
					{
						fibNode swap = x;
						x = y;
						y = swap;
					}

					Fib_Heap_Link(H, y, x);		// The tree with the lesser value becomes parent, the other becomes a child.
					A[d] = null;
					d = d + 1;
				}
				A[d] = x;						// Put either the unaltered tree or the merged tree into its proper place in the array.
				x = H.roots.removeLast();		// Keep pulling nodes until you get a null node, which indicates root list is empty.
			}
			
			H.min = null;						// Clear out min pointer and root list. Get ready to rebuild root list.
			H.roots.clear();
			
			for (int i = 0; i < A.length; i++)	// Rebuild the root list from the temporary array.
			{
				if (A[i] != null)
				{
					H.roots.add(A[i]);
					if (H.min == null || A[i].value < H.min.value)
					{
						H.min = A[i];
					}
				}
			}
		}
		
		public void Fib_Heap_Link(fibHeap H, fibNode y, fibNode x)	// Removes y from the root list of H and makes it a child of x.
		{
			
			// When this method is called, the root list is already destroyed and in the process of being consolidated.
			// So can skip the step of removing y from the root list.
			
			if (x.child == null)		// If x doesn't have children, create new child list and add y to it.
			{
				x.child = new doubleFibChain();
				y.parent = x;
				y.mark = false;
				x.child.add(y);
			}
			else						// If x already has children, add y to the child list.
			{
				y.parent = x;
				y.mark = false;
				x.child.add(y);
			}
			x.degree++;					// Increase x's degree by 1 because it has one more child than previous.
		}
		
		public void Fib_Heap_Decrease_Value(fibHeap H, fibNode x, int k)
		{
			if (k > x.value)
			{
				// Do nothing because new key is greater than current key.
			}
			else
			{
				x.value = k;		// Set x's value to the new lower value.
				x.mark = false;		// Reset x's mark to false.
				fibNode y = x.parent;
				if (y != null && x.value < y.value)		// If x is a child node, cut x from its parent and bring it to the root list.
				{
					Cut(H, x, y);
					Cascading_Cut(H, y);				// Also cascade up the tree.
				}
				if (x.value < H.min.value)
				{
					H.min = x;							// Update min pointer if necessary.
				}
			}
		}
		
		public void Cut(fibHeap H, fibNode x, fibNode y)	// "Cuts" the link between x and its parent y, making x a root.
		{
			// Remove x from the child linked list of y, decrementing y.degree.
			y.child.remove(x);
			
			y.degree--;
			
			x.parent = null;
			x.mark = false;
			
			// Add x to the root list of H and replace H.min if x is lower.
			H.roots.add(x);
			if (x.value < H.min.value)
			{
				H.min = x;
			}

		}
		
		public void Cascading_Cut(fibHeap H, fibNode y)		// Recurses up the tree until it finds either a root or unmarked node.
		{
			fibNode z = y.parent;
			if (z != null)
			{
				if (y.mark == false)			// If the node is unmarked, mark it.
				{
					y.mark = true;
				}
				else							// If the node is already marked, cut it. And test its parent.
				{
					Cut(H, y, z);
					Cascading_Cut(H, z);
				}
			}
		}
		
		public void Fib_Heap_Delete(fibHeap H, fibNode x)
		{
			Fib_Heap_Decrease_Value(H, x, Integer.MIN_VALUE);
			Extract_Min(H);
		}
		
		public boolean isEmpty(fibHeap H)
		{
			if (H.min != null)
			{
				return false;
			}
			else
			{
				return true;
			}
		}
	}
		
	class fheapScheme			// Main object for implementing the Fibonacci heap scheme for Dijkstra's algorithm. 
	{
		Graph igraph;				// The input graph
		int startindex;				// Index to the source node
		fibNode[] fibNodeCache;		// Array storage that holds pointers to each fibNode for O(1) access.
		int[] cost;					// Array to store the final cost results of dijkstra's algorithm.
		int[] previous;				// Each node's previous node for dijkstra's algorithm.
		FibonacciHeap fib;			// Control object for fibonacci heap.
		fibHeap Q;					// Priority queue using a Fibonacci Heap.
		HashSet<Integer> Qa;		// Helper data structure for quick "contains" checking.
		long runningtime;			
		
		public fheapScheme(Graph graph)
		{
			igraph = graph;
			startindex = graph.x;
			fibNodeCache = new fibNode[igraph.n];
			cost = new int[graph.n];
			previous = new int[graph.n];
			fib = new FibonacciHeap();
			Q = fib.Make_Fib_Heap();
			Qa = new HashSet<Integer>();
			runningtime = 0;
			
			if (!graph.verified)
			{
				if (!outputresultonly)
				{
					System.out.println("Graph has not been verified to be fully connected, or it has been determined it is not fully connected. Abort.");
				}
			}
			else
			{
				executeDijkstra();			// Starts the Dijkstra algorithm.
			}
		}
		
		public void executeDijkstra()
		{
			long starttime = System.currentTimeMillis();		// Starts the measurement of execution time.
			
			for (int i = 0; i < igraph.n; i++)	// Initialization step. Set all costs to infinity and previous nodes to undefined.
			{
				fibNodeCache[i] = new fibNode(i, Integer.MAX_VALUE);	// Highest integer value chosen to mean "infinity".
				fib.Fib_Heap_Insert(Q, fibNodeCache[i]);				// At start all nodes unoptimized, so all are in Q.
				Qa.add(i);
				previous[i] = -99;										// Arbitrary value < 0 chosen to mean "undefined".
			}
			
			fib.Fib_Heap_Decrease_Value(Q, fibNodeCache[startindex], 0);	// Set the distance from source to source as 0.		
			
			
			while (!fib.isEmpty(Q))	// The main loop
			{
				fibNode closest = fib.Extract_Min(Q);			// Extracts the next closest node to examine.
				Qa.remove(Integer.valueOf(closest.key));		// Removes the node's key from the helper set.
				cost[closest.key] = closest.value;				// Sets the final result of that node's cost.
				int u = closest.key;
				if (closest.value == Integer.MAX_VALUE)
				{
					break;						// Means node is unreachable, end algorithm.
				}
				
				LinkedList<Vertex> V = igraph.graphs[u].vertices;		// Get all neighbors of node u
				ListIterator<Vertex> l2a = V.listIterator();
				
				while (l2a.hasNext())		// For each neighbor v of u where v has not been removed from Q.
				{
					Vertex v = l2a.next();
					if (Qa.contains(Integer.valueOf(v.index)))	// Check if node v is still in Q via helper Qa.
					{
						int temp = closest.value + v.distance;		// Relax (u,v).
						if (temp < fibNodeCache[v.index].value)		// Update cost if necessary.
						{
							fib.Fib_Heap_Decrease_Value(Q, fibNodeCache[v.index], temp);
							previous[v.index] = u;
						}
					}
				}
				// Repeat process until there are no more nodes to examine.
			}
			
			// Calculate the runtime.
			long endtime = System.currentTimeMillis();
			runningtime = endtime - starttime;
			
			// Output the results to the console.
			if (outputfullresultpath && !outputresultonly)
			{
				for (int i = 0; i < igraph.n; i++)
				{
					System.out.print(cost[i]+"\t// cost from node "+igraph.x+" to "+i);
					int backtrack = previous[i];
					while (backtrack != -99)
					{
						System.out.print(" from "+backtrack);
						backtrack = previous[backtrack];
					}
					System.out.println("");
				}
			}
			else if (suppressresultoutput && !outputresultonly)
			{
				
			}
			else if (outputresultonly)
			{
				for (int i = 0; i < igraph.n; i++)
				{
					System.out.println(cost[i]);
				}
			}
			else
			{
				for (int i = 0; i < igraph.n; i++)
				{
					System.out.println(cost[i]+"\t// cost from node "+igraph.x+" to "+i);
				}
			}
		}
	}
	
	public static int randomi(int n1)		// Random function returns random integer from 0 to k-1.
	{
		Random generator = new Random(System.nanoTime());
		return generator.nextInt(n1);
	}
	
	public static boolean numsAndOneDecimal(String s)		// This helper method checks for proper input.
	{
		boolean decimal = false;
		for (int i = 0; i < s.length(); i++)
		{
			if (s.charAt(i) < 48 && s.charAt(i) > 57)
			{
				if (s.charAt(i) == 46 && decimal == false)
				{
					decimal = true;
				}
				else
				{
					return false;
				}
			}
		}
		return true;
	}
	
	class readFile	// This helper method reads the text file and initiates the user input modes.
	{
		public void readNodes(String filepath, int mode) throws IOException
		{
			InputStream is = null;		// Important that the filepath NOT have spaces!
			if (!outputresultonly)
			{
				System.out.println("Loading user input text file at filepath: "+filepath);
			}
			try
			{
				is = new FileInputStream(filepath);
				fastinput.init(is);
				int x = fastinput.nextInt();
				int n = fastinput.nextInt();
				int m = fastinput.nextInt();
				userInputEdges[] edges = new userInputEdges[m];
				for (int i = 0; i < m; i++)
				{
					int lv1 = fastinput.nextInt();
					int lv2 = fastinput.nextInt();
					int lc = fastinput.nextInt();
					edges[i] = new userInputEdges(lv1, lv2, lc);
				}
				if (!outputresultonly)
				{
					System.out.println("Creating new graph with user provided parameters: x = "+x+" n = "+n+" m = "+m);
					System.out.println("Edges are: ");
					for (int i = 0; i < m; i++)
					{
						System.out.println(edges[i].v1+" "+edges[i].v2+" "+edges[i].c);
					}
				}
				
				Graph s = new Graph(x, n, m, edges, 0);
				if (mode == 1)
				{
					if (!outputresultonly)
					{
						System.out.println("Starting simple scheme run.");
					}
					simpleScheme ss = new simpleScheme(s);
					if (s.verified && !outputresultonly)
					{
						System.out.println("Finished simple scheme run. Duration = "+ss.runningtime+" milliseconds.");
					}
					else if (outputresultonly)
					{
						
					}
					else
					{
						System.out.println("Aborted simple scheme run.");
					}
				}
				if (mode == 2)
				{
					if (!outputresultonly)
					{
						System.out.println("Starting f-heap scheme run.");
					}
					
					fheapScheme fh = new fheapScheme(s);
					if (s.verified && !outputresultonly)
					{
						System.out.println("Finished f-heap scheme run. Duration = "+fh.runningtime+" milliseconds.");
					}
					else if (outputresultonly)
					{
						
					}
					else
					{
						System.out.println("Aborted f-heap scheme run.");
					}
				}
			}
			catch (Exception e)
			{
				System.out.println("Unable to access the file (does not exist at file path or cannot read). Please try again.");
			}
			finally
			{
				if (is != null)
				{
					is.close();
				}
			}
		}
	}
	
	public static class fastinput		// This helper method reads the text file.
	{
		static BufferedReader bf;
		static StringTokenizer st;
		
		static void init(InputStream fr)
		{
			bf = new BufferedReader(new InputStreamReader(fr));
			st = new StringTokenizer("");
		}
		
		static String nextline() throws IOException
		{
			return bf.readLine();
		}
		
		static String next() throws IOException
		{
			while (!st.hasMoreTokens())
			{
				st = new StringTokenizer(bf.readLine());
			}
			return st.nextToken();
		}
		
		static int nextInt() throws IOException
		{
			return Integer.parseInt(next());
		}
	}
	
	public static void main(String[] args) throws IOException		// Main method.
	{
		if (args.length == 0)
		{
			System.out.println("Error: please input an argument.");
		}
		else if (args.length == 2 && args[0].equals("-s"))
		{
			dijkstra c = new dijkstra();
			readFile r = c.new readFile();
			r.readNodes(args[1], 1);
		}
		else if (args.length == 2 && args[0].equals("-f"))
		{
			dijkstra c = new dijkstra();
			readFile r = c.new readFile();
			r.readNodes(args[1], 2);
		}
		else if (args.length == 4 && args[0].equals("-r"))
		{
			boolean validinput = true;
			for (int i = 1; i < 4; i++)
			{
				if (!numsAndOneDecimal(args[i]))
				{
					validinput = false;
					break;
				}
			}
			
			if (!validinput)
			{
				System.out.println("Invalid argument parameters. Please try again.");
			}
			else if (args[1].equals("1000") && args[2].equals("0.1"))
			{
				System.out.println("Test case of n = 1000 and d = 0.1% is to be skipped!");
			}
			else
			{
				int z = 1;
				long totalrunningsimple = 0;
				long totalrunningfib = 0;
				int simplebetter = 0;
				int fheapbetter = 0;
				int ties = 0;
				long shortsimple = Long.MAX_VALUE;
				long longsimple = Long.MIN_VALUE;
				long shortfheap = Long.MAX_VALUE;
				long longfheap = Long.MIN_VALUE;
				long[] sstimes = new long[10];
				long[] fhtimes = new long[10];
//				for (;z <= 10; z++)
//				{
					dijkstra b = new dijkstra();
//					System.out.println("");
//					System.out.println("Run # "+z+" of 10:");
					Graph a = b.new Graph(Integer.valueOf(args[1]), Double.valueOf(args[2]), Integer.valueOf(args[3]), 1);
//					System.out.println("Starting simple scheme run.");
					simpleScheme ss = b.new simpleScheme(a);
//					System.out.println("Finished simple scheme run.");
					System.out.println("");
//					System.out.println("Starting Fibonacci scheme run.");
					fheapScheme fs = b.new fheapScheme(a);
//					System.out.println("Finished Fibonacci scheme run.");
//					System.out.println("");
//					System.out.println("Simple scheme run duration = "+ss.runningtime+" milliseconds.");
//					System.out.println("F-heap scheme run duration = "+fs.runningtime+" milliseconds.");
					sstimes[z-1] = ss.runningtime;
					fhtimes[z-1] = fs.runningtime;
					totalrunningsimple = totalrunningsimple + ss.runningtime;
					totalrunningfib = totalrunningfib + fs.runningtime;
					
					if (ss.runningtime < fs.runningtime)
					{
						simplebetter++;
					}
					else if (ss.runningtime > fs.runningtime)
					{
						fheapbetter++;
					}
					else
					{
						ties++;
					}
					if (ss.runningtime < shortsimple)
					{
						shortsimple = ss.runningtime;
					}
					if (ss.runningtime > longsimple)
					{
						longsimple = ss.runningtime;
					}
					if (fs.runningtime < shortfheap)
					{
						shortfheap = fs.runningtime;
					}
					if (fs.runningtime > longfheap)
					{
						longfheap = fs.runningtime;
					}
					
					System.out.println("");
//				}
//				System.out.println("Finished all 10 runs of simple scheme and f-heap scheme.");
//				System.out.println("");
//				System.out.println("Simple scheme average duration (total/10) = "+(totalrunningsimple/10)+" milliseconds.");
//				System.out.println("F-heap scheme average duration (total/10) = "+(totalrunningfib/10)+" milliseconds.");
//				System.out.println("");
//				System.out.println("Simple scheme was faster "+simplebetter+" times. F-heap scheme was faster "+fheapbetter+" times with "+ties+" ties.");
//				System.out.println("Improved average results (lowest and highest times omitted with average of remaining times):");
				long improvedsimple = totalrunningsimple - shortsimple - longsimple;
				long improvedfheap = totalrunningfib - shortfheap - longfheap;
//				System.out.println("Simple scheme = "+(improvedsimple/8)+" milliseconds.  F-heap scheme = "+(improvedfheap/8)+" milliseconds.");
				Arrays.sort(sstimes);
				Arrays.sort(fhtimes);
				long simplemed = (sstimes[4]+sstimes[5])/2;
				long fheapmed = (fhtimes[4]+fhtimes[5])/2;
//				System.out.println("Simple scheme median time = "+simplemed+" milliseconds.  F-heap scheme median time = "+fheapmed+" milliseconds.");
//				System.out.println("");
			}
		}
		else
		{
			System.out.println("Invalid argument parameters. Please try again.");
		}
	}
}
