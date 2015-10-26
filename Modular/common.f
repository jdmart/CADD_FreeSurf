	Module mod_Common
	TYPE region
	double precision xmin, xmax, ymin, ymax, damped_width
	END TYPE region
	

	TYPE MD_Thermostat
	Character(len=16) :: Type
	Character (len=16) :: Damping_Mode
	end TYPE MD_Thermostat			
	end module mod_common
