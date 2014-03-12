%global pMATLAB
pMatlab_Init;

script_file_head = '.pRUN_Parallel_Stub_';
script_file_tail = '_temp';
script_file = dir([script_file_head,'*',script_file_tail]);
m_file =  regexp(script_file.name, [script_file_head,'(.*)',script_file_tail], 'tokens');

clear('script_file_head','script_file','script_file_tail');

m_file = char(m_file{1});

eval(m_file);

pMatlab_Finalize;
pMATLAB.my_rank = 0;
pMATLAB.comm_size = 1;


