nb = 6;
len_motif = 5; % length of motif to find
x = '%05.f';

bpos_now = bpos(:,nb);
[bland, btoff] = findchanges(bflying(:,nb),1,1);
btoff = find(btoff); btoff = btoff(2:end);
bland = find(bland); bland = bland(1:end-1);


bpos_now_corrected = bpos_now;
pos_values = unique(bpos_now);
for i = 1:length(btoff)
    [u,~,uidx] = unique(bpos_now(bland(i):btoff(i)));
    counts = accumarray(uidx,1);
    bpos_now_corrected(bland(i):btoff(i)) = u(counts == max(counts));
end


landings = [bland,btoff,bpos_now(bland),bpos_now(btoff),bpos_now_corrected(bland),bpos_now_corrected(btoff)];


pos_seq = bpos_now_corrected(bland);
u_pos = unique(pos_seq); % possible positions


seq_list = str2num(dec2base(0:length(u_pos)^len_motif-1,length(u_pos)));
seq_list = num2str(seq_list,x);
pot_seq = zeros(size(seq_list)); % potential sequence
for i = 1:size(pot_seq,1)
    for j = 1:size(pot_seq,2)
        pot_seq(i,j) = u_pos(str2num(seq_list(i,j))+1);
    end
end

seq_count = zeros(size(pot_seq,1),1); % count of sequence motif
for i = 1:length(pos_seq)-len_motif+1
    idx = find(ismember(pot_seq,pos_seq(i:i+len_motif-1)','rows'));
    seq_count(idx) = seq_count(idx) + 1;
end


