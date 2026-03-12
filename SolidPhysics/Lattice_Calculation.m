function Crystal = Lattice_Calculation(Conventional_size, Conventional_angle, lattice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 单位(unit)：(size)/Å (angle)/rad
% input:  [Conventional_size] [Conventional_angle]
% output: cell{CrystalSystem_type, lattice_type, atom_position, Primitive_size, Primitive_angle}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    Conventional_size (:,1) double
    Conventional_angle (:,1) double
    lattice (1,:) double
end
tol = 1e-9;
[a_Con, b_Con, c_Con] = deal(Conventional_size(1), Conventional_size(2), Conventional_size(3));
[alpha_Con, beta_Con, gamma_Con] = deal(Conventional_angle(1), Conventional_angle(2), Conventional_angle(3));
[P, I, C, F] = deal(lattice(1), lattice(2), lattice(3), lattice(4));
CrystalSystem_type = [[abs(a_Con - b_Con) < tol, abs(a_Con - c_Con) < tol, abs(b_Con - c_Con) < tol], ...
    [abs(alpha_Con - beta_Con) < tol, abs(alpha_Con - gamma_Con) < tol, abs(beta_Con - gamma_Con) < tol], ...
    [abs(alpha_Con - pi/2) < tol, beta_Con < 2*pi/3, abs(gamma_Con - 2*pi/3) < tol]];

if all(CrystalSystem_type == [[1 1 1], [1 1 1], [1 1 0]])
    CrystalSystem_type = "Cubic";
    if P > 0
        lattice_type = "cP";
        atom_position = [0; 0; 0];
        a_Pri = a_Con; b_Pri = b_Con; c_Pri = c_Con;
        alpha_Pri = alpha_Con; beta_Pri = beta_Con; gamma_Pri = gamma_Con;
        a1_Pri = [a_Con; 0; 0]; a2_Pri = [0; b_Con; 0]; a3_Pri = [0; 0; c_Con];
    elseif I > 0
        lattice_type = "cI";
        atom_position = [[0; 0; 0], [0.5; 0.5; 0.5]];
        a_Pri = sqrt(3)*0.5*a_Con; b_Pri = sqrt(3)*0.5*b_Con; c_Pri = sqrt(3)*0.5*c_Con;
        alpha_Pri = acos(-1/3); beta_Pri = acos(-1/3); gamma_Pri = acos(-1/3);
        a1_Pri = 0.5*[a_Con; b_Con; -c_Con]; a2_Pri = 0.5*[-a_Con; b_Con; c_Con]; a3_Pri = 0.5*[a_Con; -b_Con; c_Con];
    elseif C > 0
        lattice_type = "cC";
        error(["非标准类型，该点阵类型不存在",lattice_type])
    elseif F > 0
        lattice_type = "cF";
        atom_position = [[0; 0; 0], [0.5; 0.5; 0], [0.5; 0; 0.5], [0; 0.5; 0.5]];
        a_Pri = sqrt(2)*0.5*a_Con; b_Pri = sqrt(2)*0.5*b_Con; c_Pri = sqrt(2)*0.5*c_Con;
        alpha_Pri = acos(1/2); beta_Pri = acos(1/2); gamma_Pri = acos(1/2);
        a1_Pri = 0.5*[0; b_Con; c_Con]; a2_Pri = 0.5*[a_Con; 0; c_Con]; a3_Pri = 0.5*[a_Con; b_Con; 0];
    end
elseif all(CrystalSystem_type == [[1 0 0], [1 0 0], [1 1 1]])
    CrystalSystem_type = "Hexagonal";
    if P > 0
        lattice_type = "hP";
        atom_position = [[0; 0; 0], [2/3; 1/3; 0]];
        a_Pri = a_Con; b_Pri = b_Con; c_Pri = c_Con;
        alpha_Pri = alpha_Con; beta_Pri = beta_Con; gamma_Pri = gamma_Con;
        a1_Pri = [a_Con; 0; 0]; a2_Pri = [a_Con; cos(gamma_Con)*b_Con; 0]; a3_Pri = 0.5*[0; 0; c_Con];
    elseif I > 0
        lattice_type = "hI";
        error(["非标准类型，该点阵类型不存在",lattice_type])
    elseif C > 0
        lattice_type = "hC";
        error(["非标准类型，该点阵类型不存在",lattice_type])
    elseif F > 0
        lattice_type = "hF";
        error(["非标准类型，该点阵类型不存在",lattice_type])
    end
elseif all(CrystalSystem_type == [[1 1 1], [1 1 1], [0 1 0]])
    CrystalSystem_type = "Trigonal";
    if P > 0
        lattice_type = "hR";
        atom_position = [0; 0; 0];
        a_Pri = a_Con; b_Pri = b_Con; c_Pri = c_Con;
        alpha_Pri = alpha_Con; beta_Pri = beta_Con; gamma_Pri = gamma_Con;
        a1_Pri = [a_Con; 0; 0]; a2_Pri = 0.5*[a_Con; 0; c_Con]; a3_Pri = 0.5*[a_Con; b_Con; 0];
    else
        error("该点阵类型不存在")
    end
elseif all(CrystalSystem_type == [[1 0 0], [1 1 1], [1 1 0]])
    CrystalSystem_type = "Tetragonal";
    if P > 0
        lattice_type = "tP";
        atom_position = [0; 0; 0];
        a_Pri = a_Con; b_Pri = b_Con; c_Pri = c_Con;
        alpha_Pri = alpha_Con; beta_Pri = beta_Con; gamma_Pri = gamma_Con;
        a1_Pri = [a_Con; 0; 0]; a2_Pri = [0; b_Con; 0]; a3_Pri = [0; 0; c_Con];
    elseif I > 0
        lattice_type = "tI";
        atom_position = [[0; 0; 0], [0.5; 0.5; 0.5]];
        a_Pri = a_Con; b_Pri = sqrt(2*a_Con^2+c_Con^2)*0.5; c_Pri = c_Con;
        alpha_Pri = acos(c_Con/sqrt(2*a_Con^2+c_Con^2)); beta_Pri = beta_Con; gamma_Pri = acos(a_Con/sqrt(2*a_Con^2+c_Con^2));
        a1_Pri = [a_Con; 0; 0]; a2_Pri = [a_Con; b_Con; c_Con]; a3_Pri = [0; 0; c_Con];
    elseif C > 0
        lattice_type = "tC";
        error(["非标准类型，该点阵类型不存在",lattice_type])
    elseif F > 0
        lattice_type = "tF";
        error(["非标准类型，该点阵类型不存在",lattice_type])
    end
elseif all(CrystalSystem_type == [[0 0 0], [1 1 1], [1 1 0]])
    CrystalSystem_type = "Orthorhombic";
    if P > 0
        lattice_type = "oP";
        atom_position = [0; 0; 0];
        a_Pri = a_Con; b_Pri = b_Con; c_Pri = c_Con;
        alpha_Pri = alpha_Con; beta_Pri = beta_Con; gamma_Pri = gamma_Con;
        a1_Pri = [a_Con; 0; 0]; a2_Pri = [0; b_Con; 0]; a3_Pri = [0; 0; c_Con];
    elseif I > 0
        lattice_type = "oI";
        atom_position = [[0; 0; 0], [0.5; 0.5; 0.5]];
        a_Pri = a_Con; b_Pri = sqrt(a_Con^2+b_Con^2+c_Con^2)*0.5; c_Pri = c_Con;
        alpha_Pri = acos(c_Con/sqrt(a_Con^2+b_Con^2+c_Con^2)); beta_Pri = beta_Con; gamma_Pri = acos(a_Con/sqrt(a_Con^2+b_Con^2+c_Con^2));        a1_Pri = [a_Con; 0; 0]; a2_Pri = [a_Con; b_Con; c_Con]; a3_Pri = [0; 0; c_Con];
        a1_Pri = [a_Con; 0; 0]; a2_Pri = [a_Con; b_Con; c_Con]; a3_Pri = [0; 0; c_Con];
    elseif C > 0
        lattice_type = "oC";
        atom_position = [[0; 0; 0], [0.5; 0.5; 0]];
        a_Pri = a_Con; b_Pri = sqrt(a_Con^2+b_Con^2)*0.5; c_Pri = c_Con;
        alpha_Pri = alpha_Con; beta_Pri = beta_Con; gamma_Pri = acos(a_Con/sqrt(a_Con^2+b_Con^2));
        a1_Pri = [a_Con; 0; 0]; a2_Pri = [a_Con; b_Con; 0]; a3_Pri = [0; 0; c_Con];
    elseif F > 0
        lattice_type = "oF";
        atom_position = [[0; 0; 0], [0.5; 0.5; 0], [0.5; 0; 0.5], [0; 0.5; 0.5]];
        a_Pri = sqrt(a_Con^2+c_Con^2)*0.5; b_Pri = sqrt(a_Con^2+b_Con^2)*0.5; c_Pri = sqrt(b_Con^2+c_Con^2)*0.5;
        alpha_Pri = acos(b_Con^2/sqrt((a_Con^2+b_Con^2)*(b_Con^2+c_Con^2))); beta_Pri = acos(c_Con^2/sqrt((b_Con^2+c_Con^2)*(a_Con^2+c_Con^2))); gamma_Pri = acos(a_Con^2/sqrt((a_Con^2+c_Con^2)*(a_Con^2+b_Con^2)));
        a1_Pri = 0.5*[0; b_Con; c_Con]; a2_Pri = 0.5*[a_Con; 0; c_Con]; a3_Pri = 0.5*[a_Con; b_Con; 0];
    end
elseif all(CrystalSystem_type == [[0 0 0], [0 1 0], [1 0 0]]) || all(CrystalSystem_type == [[0; 0; 0], [0 1 0], [1 1 0]])
    CrystalSystem_type = "Monoclinic";
    if P > 0
        lattice_type = "mP";
        atom_position = [0; 0; 0];
        a_Pri = a_Con; b_Pri = b_Con; c_Pri = c_Con;
        alpha_Pri = alpha_Con; beta_Pri = beta_Con; gamma_Pri = gamma_Con;
        a1_Pri = [a_Con; 0; 0]; a2_Pri = [0; cos(gamma_Pri)*b_Con; 0]; a3_Pri = [0; 0; c_Con];
    elseif I > 0
        lattice_type = "mI";
        atom_position = [[0; 0; 0], [0.5; 0.5; 0.5]];
        a_Pri = a_Con; b_Pri = sqrt(a_Con^2+b_Con^2+c_Con^2)*0.5; c_Pri = c_Con;
        alpha_Pri = acos(c_Con/sqrt(a_Con^2+b_Con^2+c_Con^2)); beta_Pri = beta_Con; gamma_Pri = acos(a_Con/sqrt(a_Con^2+b_Con^2+c_Con^2));
        a1_Pri = [a_Con; 0; 0]; a2_Pri = [a_Con; cos(gamma_Pri)*b_Con; c_Con]; a3_Pri = [0; 0; c_Con];
    elseif C > 0
        lattice_type = "mC";
        error(["非标准类型，该点阵类型不存在",lattice_type])
    elseif F > 0
        lattice_type = "mF";
        error(["非标准类型，该点阵类型不存在",lattice_type])
    end
elseif all(CrystalSystem_type == [[0 0 0], [0 0 0], [0 0 0]]) || all(CrystalSystem_type == [[0 0 0], [0 0 0], [0 1 0]]) || all(CrystalSystem_type == [[0 0 0], [0 0 0], [1 0 0]]) || all(CrystalSystem_type == [[0 0 0], [0 0 0], [1 1 0]])
    CrystalSystem_type = "Triclinic";
    if P > 0
        lattice_type = "aP";
        atom_position = [0; 0; 0];
        a_Pri = a_Con; b_Pri = b_Con; c_Pri = c_Con;
        alpha_Pri = alpha_Con; beta_Pri = beta_Con; gamma_Pri = gamma_Con;
    elseif I > 0
        lattice_type = "aI";
        error(["非标准类型，该点阵类型不存在",lattice_type])
    elseif C > 0
        lattice_type = "aC";
        error(["非标准类型，该点阵类型不存在",lattice_type])
    elseif F > 0
        lattice_type = "aF";
        error(["非标准类型，该点阵类型不存在",lattice_type])
    end
else
    error("非标准晶格常数输入，布拉菲点阵不存在，请检查")
end
Primitive_size = [a_Pri; b_Pri; c_Pri];
Primitive_angle = [alpha_Pri; beta_Pri; gamma_Pri];
Primitive_cell = [a1_Pri, a2_Pri, a3_Pri];
Crystal = {CrystalSystem_type, lattice_type, atom_position, Primitive_size, Primitive_angle, Primitive_cell};
end

