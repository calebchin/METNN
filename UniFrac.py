def weighted_unifrac(P, Q, L):
    """
    Compute the Weighted UniFrac distance between two samples P and Q.

    Parameters:
    - P: Empirical distribution (relative abundance) of sample P
    - Q: Empirical distribution (relative abundance) of sample Q
    - L: List of edge lengths associated with the OTUs in the tree

    Returns:
    - Weighted UniFrac distance
    """
    # Initialize the numerator of the Weighted UniFrac distance
    numerator = 0

    #loop through all OTU, calc the abs difference between the relative abudance
    #multiply that by edge length

    for i in range(len(P)):
        numerator += abs(P[i] - Q[i]) * L[i]

    #total edgelength
    total_edge_length = sum(L)

    #  Weighted UniFrac distance
    weighted_unifrac_distance = numerator / total_edge_length

    return weighted_unifrac_distance