"""Small ctypes adapter for the production C ABI and compatibility gates."""
import ctypes as C
from pathlib import Path
class Engine:
    def __init__(self,library,system):
        self.L=C.CDLL(str(library));L=self.L;h=C.c_void_p;dp=C.POINTER(C.c_double);ip=C.POINTER(C.c_int)
        signatures={
          'set_log_directory':(None,[C.c_char_p]),'set_log_level':(None,[C.c_int]),
          'create_from_lst':(h,[C.c_char_p,C.c_char_p,C.c_int]),'destroy':(None,[h]),
          'num_elements':(C.c_int,[h]),'num_species':(C.c_int,[h]),'num_phases':(C.c_int,[h]),
          'element_name':(C.c_char_p,[h,C.c_int]),'species_name':(C.c_char_p,[h,C.c_int]),
          'phase_name':(C.c_char_p,[h,C.c_int]),'element_index':(C.c_int,[h,C.c_char_p]),
          'species_index':(C.c_int,[h,C.c_char_p]),'phase_index':(C.c_int,[h,C.c_char_p]),
          'species_phase':(C.c_int,[h,C.c_int]),'species_is_gas':(C.c_int,[h,C.c_int]),
          'species_molar_mass':(C.c_double,[h,C.c_int]),'stoich':(C.c_double,[h,C.c_int,C.c_int]),
          'TP_range':(None,[h,dp,dp,dp,dp]),'state_size':(C.c_int,[h]),
          'set_species_bounds':(None,[h,C.c_int,C.c_double,C.c_double]),
          'equilibrate':(C.c_int,[h,C.c_double,C.c_double,dp,C.c_int,dp,ip]),
          'species_amounts':(None,[h,dp]),'gas_partial_pressures':(None,[h,dp]),
          'log10_activities':(None,[h,dp]),'phase_amounts':(None,[h,dp]),
          'phase_log10_saturation':(None,[h,dp]),'last_iterations':(C.c_int,[h]),
          'last_seconds':(C.c_double,[h]),'last_error':(C.c_char_p,[h]),
          'suppress_condensed':(None,[h,C.c_int]),'species_g0':(C.c_int,[h,C.c_double,C.c_double,dp]),
          'element_potentials':(None,[h,dp]),'gas_mole_fractions':(None,[h,dp]),
          'pair_peq':(C.c_int,[h,C.c_int,ip,C.c_int,C.c_int,dp,ip]),
          'balance_error':(C.c_double,[h,dp,C.c_double]),
          'set_ipm_controls':(None,[h,C.c_double,C.c_double,C.c_int]),
          'get_ipm_controls':(None,[h,dp,dp,ip]),'set_warm_iteration_limit':(None,[h,C.c_int])}
        for n,(ret,args) in signatures.items():
            if not hasattr(L,'gemsb_'+n):continue
            f=getattr(L,'gemsb_'+n);f.restype=ret;f.argtypes=args
        L.gemsb_set_log_directory(b'/tmp');L.gemsb_set_log_level(6)
        error=C.create_string_buffer(2048);self.h=L.gemsb_create_from_lst(str(system).encode(),error,len(error))
        assert self.h,error.value
        self.elements=[L.gemsb_element_name(self.h,i).decode() for i in range(L.gemsb_num_elements(self.h))]
        self.species=[L.gemsb_species_name(self.h,j).decode() for j in range(L.gemsb_num_species(self.h))]
        self.is_gas=[L.gemsb_species_is_gas(self.h,j) for j in range(len(self.species))]
        self.state=(C.c_double*L.gemsb_state_size(self.h))();self.valid=C.c_int(0)
    def close(self):
        if self.h:self.L.gemsb_destroy(self.h);self.h=None
    def suppress(self):
        for j,g in enumerate(self.is_gas):
            if not g:self.L.gemsb_set_species_bounds(self.h,j,0,0)
    def equilibrate(self,T,P,b,mode=2):
        self.bulk=(C.c_double*len(self.elements))(*(b.get(e,0) for e in self.elements))
        return self.L.gemsb_equilibrate(self.h,T,P,self.bulk,mode,self.state,C.byref(self.valid))
    def values(self,name,size=None):
        v=(C.c_double*(size or len(self.species)))();getattr(self.L,'gemsb_'+name)(self.h,v);return list(v)
    def g0(self,T,P):
        v=(C.c_double*len(self.species))();assert self.L.gemsb_species_g0(self.h,T,P,v)==0;return list(v)
