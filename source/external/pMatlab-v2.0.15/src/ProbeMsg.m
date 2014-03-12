function [message_rank, numeric_tag, string_tag] = ProbeMsg(source, tag)
%
global pMATLAB

[message_rank, numeric_tag, string_tag] = MPI_Probe(source, tag, pMATLAB.comm);

