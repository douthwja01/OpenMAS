

trajectories = rand([36,500]);


spmd
  D = codistributed(trajectories)
end