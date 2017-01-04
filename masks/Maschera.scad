	angoloTast = -5;
	angoloSpall = -5;
	angoloClavicola = -15;
	angoloFresa = -5;
	fermo = 35+4;
	pollice = 10.6;
	fermo2 = fermo + pollice;
	riccio = 14.5;
	fresa = 9;
	stretto = 2.3;
	largo = 4.2;
	lungo = 27;
	cassaInit = 4;
	cassaLungo = 21;
	cassaRaggio = 5;
	cassaSpessore=0.4;
	asta = cassaLungo+cassaInit;
	
	module manico(){
		difference(){
		union() {
		
		//primo tratto
			translate([asta,0,0])cube([fermo-asta, 2, 2]);
			translate([0,1,1])rotate([0,90,0])cylinder(h=asta,d=0.8,$fn=20);
		//secondo tratto
			translate([fermo2,0,0])	rotate([0,angoloFresa,0])cube([riccio,2,4]);
		//fermi
			translate([fermo-2,0,0])cube([2,2,4]);
			translate([fermo2,0,0])	cube([2,2,4]);
			translate([fermo-4,0,1.9])rotate([0,45,0])cube([2,2,3]);
			translate([fermo2+2.5,0,0.5])rotate([0,-45,0])cube([2,2,3]);
		//pollice
			translate([fermo,1,1])rotate([0,90,0])cylinder(h=pollice,d=2,$fn=20);
		}
		union(){
			translate([fermo,0,-1])rotate([0,angoloTast,0])	cube([pollice,2,0.5]);
			translate([fermo2+2,0.5,-.5])rotate([0,angoloTast,0])cube([fresa,1,3]);
		//piroli
			translate([fermo2+3,2.5,1.6])rotate([90,0,0])cylinder(r1=0.6,r2=0.8,h=3,$fn=20);
			translate([fermo2+5,-0.5,1.8])rotate([-90,0,0])cylinder(r1=0.6,r2=0.8,h=3,$fn=20);
			translate([fermo2+7.5,2.5,2])rotate([90,0,0])cylinder(r1=0.6,r2=0.8,h=3,$fn=20);
			translate([fermo2+9.5,-0.5,2.2])rotate([-90,0,0])cylinder(r1=0.6,r2=0.8,h=3,$fn=20);	
		//becco
			translate([fermo2+6,-1,19.5])rotate([-90,0,0])cylinder(r=16,h=4,$fn=40);	
		
		}
		}
	}
	module tastiera(){
		difference(){
		intersection(){
		translate([0,0,0])rotate([0,-90,0])cylinder(r1=stretto,r2=largo,h=lungo,$fn=20);
		translate([0,0,4])rotate([0,-90,0])cylinder(r1=stretto,r2=largo,h=lungo,$fn=20);
		}
		union(){
		//translate([-lungo/2,0,4])
			//cube([lungo+1,2*largo,largo],center=true);
		translate([0,0,5])rotate([0,-90,0])cylinder(r1=stretto,r2=largo,h=lungo+1,$fn=20);
	
		}
		}
	}
	module spalliera(){
		difference(){
		union(){
		intersection(){
			cube([20,2,12]);
			translate([2,3,-20])rotate([90,0,0])cylinder(r=26,h=4,$fn=40);	
		}
		intersection(){
			cube([20,2,12]);
			translate([0,3,31])rotate([90,0,0])cylinder(r=26,h=4,$fn=40);	
		}
		rotate([-angoloClavicola,0,0])difference(){
		translate([0,0.5,-0.25])cube([10,6,1]);
		translate([5,8.5,-1.25])cylinder(r=6,h=4,$fn=40);	
		}
		}
		union(){
			translate([8,3,36])rotate([90,0,0])cylinder(r=26,h=4,$fn=40);	
		}
		}
	}
	module cassa(){
		difference(){
			union(){
			cylinder(r=cassaRaggio,h=cassaLungo,$fn=40);
			}
			union(){
			translate([0,0,-1])cylinder(r=cassaRaggio-cassaSpessore,h=cassaLungo+2,$fn=40);
			}
		}
	}
	module arco(){
		lato = 6;
		lato2 = 8;
		intersection(){//col cubo
		difference(){
				translate([0,0,-1])cylinder(r=lato,h=4,$fn=40);
				translate([-1.5,0,-1])cylinder(r=lato,h=4,$fn=40);
			}
			translate([-5,0,0])cube([13,lato,2]);
		}
	}

	module sagoma(){
		lato = 6;
		lato2 = 8;
		union(){
			arco();
			translate([-lato-.5,lato-.75,0])rotate([0,0,270])arco();
		}
	}
	
	color([1.0,0.7,0.1,1.0])translate([-40,0,0])manico();
	color([0.8,0.1,0.7,1.0]){
		translate([9.8,1,-1.5])rotate([0,angoloTast,0])tastiera();
	}
	color([0.2,0.1,0.7,1.0]){
		translate([-40,4,-2])difference(){rotate([angoloClavicola,0,-90])spalliera();manico();}
	}
	color([0.8,0.8,0.1,1.0]){
		translate([-fermo+cassaInit,1,4])rotate([0,90,0])cassa();
	}
	color([0.2,0.8,0.4,1.0]){
		translate([-7,2,0])rotate([0,0,0])sagoma();
	}
