function DFS(i,nnmax,nn)
    global n_marker nnode_in_cluster
    if (n_marker(i)==0)
	n_marker(i)=1;
	nnode_in_cluster=nnode_in_cluster+1;
        for j=1:nnmax(i)
	    if (n_marker(nn(j,i))==0)
	        DFS(nn(j,i),nnmax,nn);
	    end
	end
    end
end
