function prev_start_time = prevBlcokStartTime(seq)
total_dur = sum(seq.blockDurations);
prev_blk_ind = length(seq.blockEvents);
prev_blk = seq.getBlock(prev_blk_ind);
prev_blk_dur = prev_blk.blockDuration;
prev_start_time = total_dur - prev_blk_dur;

toktoktok_time = 1030e-6;
prev_start_time = prev_start_time + toktoktok_time;
end