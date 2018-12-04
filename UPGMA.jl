

function DistanceMatrix()
    f= open("input.txt")
    lines=readstring(f)

    c=parse(Float64,filter(x -> !isspace(x), lines))

    n = c[1]
    c = c[2:length(c)]
    dm = reshape(c,(n,n))
    return dm
end

function upgma(dm)
    n = length(dm[1,:])
    # each node is placed in a cluster
    # heights of leaf nodes are all set to zero
    # newick string starts off as a list of nodes
    clusters = Array{Int64}[]
    newick = []
    heights = []
    nodes = []
    for i in 1:n
        push!(clusters,[i])
        newick = vcat(newick,"$i")
        nodes = vcat(nodes,"$(i-1)")
        heights = vcat(heights,0)
    end

    next = n+1

    while n > 1
        # add 2 * max to the diagonal zeros before finding the indices and value of the min distance
        (max,ind)= findmax(dm)
        dme = dm + eye(n,n)*max*2
        (min,ind) = findmin(dme)
        # store the indicies of the min as row and col
        row = ((ind-1)%n)+1
        col = div(ind-1,n)+1

        ncr = length(clusters[row])
        ncc = length(clusters[col])
        # get distance to new cluster formula
        # append new row and new column to distance matrix
        newRow = ( ncr * dm[row,:] +
        ncc * dm[col,:] ) / (ncr + ncc)
        dm = vcat(dm,newRow')
        newCol = ( ncr * dm[:,row] +
        ncc * dm[:,col] ) / (ncr + ncc)
        dm = hcat(dm,newCol)
        # set the diag element of new row and new col to zero
        dm[n+1,n+1] = 0.0

        # append the new cluster to cluster list
        push!(clusters,vcat(clusters[row],
        clusters[col]))
        # compute height for the new cluster
        # and generate the Newick representation
        # for the new cluster
        h=min/2
        hr = ":$(@sprintf("%.3f",
        (h-heights[row])))"
        hc = ":$(@sprintf("%.3f",
        (h-heights[col])))"
        newNode = "("*newick[row]*hr*","*newick[col]*hc*")$next"
        # append the new newick rep,
        # the new height and the new node name
        # to the appropriate lists
        newick = vcat(newick,newNode)
        heights = vcat(heights,h)
        nodes = vcat(nodes,next-1)

        # make use of daleteat to remove
        # row and col items from each list
        if (row < col)
            deleteat!(clusters,[row,col])
            deleteat!(newick,[row,col])
            deleteat!(heights,[row,col])
            deleteat!(nodes,[row,col])
        else
            deleteat!(clusters,[col,row])
            deleteat!(newick,[col,row])
            deleteat!(heights,[col,row])
            deleteat!(nodes,[col,row])
        end

        # finally remove row and col
        # rows and columns
        # from the distance matrix
        dm = dm[setdiff(1:n+1,[row,col]),:]
        dm = dm[:,setdiff(1:n+1,[row,col])]
        # by now n should drop by one.
        n = length(dm[1,:])
        # increment the next node label
        next=next+1
    end
    # after the while loop
    # there should be one string in the
    # newick list representing the whole tree
    newick[1]
end
