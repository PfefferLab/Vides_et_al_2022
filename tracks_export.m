function T = tracks_export(btch, outfilename)
%eval('btch = currentBatch')
tracks = getfield(getfield(btch, 'results'),'tracks');

spots_total = 0;
for i=1:length(tracks)
	n = size(tracks{i});
	spots_total= spots_total + n(1);	
end 

out = zeros(spots_total, 7);  
spots_total = 0;
for i=1:length(tracks)
	%disp(i)
	n = size(tracks{i});
	n = n(1);
	if n>0
    	spt = tracks{i};
        out((spots_total+1):(spots_total+n),1) = i; %track number
        out((spots_total+1):(spots_total+n),2:7) = spt(:,1:6);
    	spots_total = spots_total+n;
	end
end 
T = array2table(out, 'VariableNames', {'TrackID','Frame','X','Y','A','BG','Sigma'}) 
T2 = groupsummary(T, ["Frame"],"sum","A")

T = join(T, T2)
T = renamevars(T, ["sum_A","GroupCount"], ["sum_A_allspotsInFrame","Num_spotsInFrame"])

writetable(T,outfilename)
end