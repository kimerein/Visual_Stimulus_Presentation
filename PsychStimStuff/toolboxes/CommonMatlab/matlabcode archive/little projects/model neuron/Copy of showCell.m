% draw synapses for list of cells
function showCell(cells,xLocs,yLocs,M)
	scale=100;

	numcells=length(cells);

	outMap=repmat(zeros(scale),numcells,1);
	inMap=repmat(zeros(scale),numcells,1);
numcells = 2 %BA hack
	for j=1:numcells-1

		for i=1:size(M,1)
			outMap(ceil(scale*xLocs(i))+(j-1)*scale,ceil(scale*yLocs(i)))=M(i,cells(j));
end
		for i=1:size(M,1)
			inMap(ceil(scale*xLocs(i))+(j-1)*scale,ceil(scale*yLocs(i)))=M(cells(j),i);
end
    	end
	X=[inMap outMap];
	X(1,1)=max(max(abs(X)))+1;
	X(1,2)=-1*X(1,1)-1;
	X(X==0)=X(1,1);
	X(scale*[1 2 3],:)=X(1,2);
	X(:,scale)=X(1,2);
	imagesc(X);
    title('Presynaptic    Postsynaptic')
	colormap('jet');
	newC=colormap;
	newC(1,:)=[0 0 0];
	newC(size(colormap,1),:)=[1 1 1];
	colormap(newC);
end
