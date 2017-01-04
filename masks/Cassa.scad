		lungo = 210;
		raggio = 50;
		spessore=4;
		fattore = 1.1;
		eccent = 5;
	asta = 300;
	module cassa(){
		difference(){
			union(){
			cylinder(r1=raggio,r2=fattore*raggio,h=lungo,$fn=40);
			}
			union(){
			translate([0,0,-1])cylinder(r1=raggio-spessore,r2=fattore*(raggio-spessore),h=lungo+2,$fn=40);
			}
		}
	}

translate([0,0,-lungo/2])cassa();
translate([0,raggio-spessore-4,-asta/2])rotate([0,0,0])cylinder(h=asta,d=8,$fn=20);
