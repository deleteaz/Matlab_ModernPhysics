function info = Materials_Database(material)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 单位(unit)：(size)/Å (angle)/rad
% input:  "material name"
% output: cell{atom_size, atom_angle, atom_number}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    material string
end
switch material
    case "Po"
        atom_size = 3.35 * ones(3,1);
        atom_angle =  90 * ones(3,1);
        atom_number = 84;
        [P, I, C, F] = deal(1, 0, 0, 0);
    case "Li"
        atom_size = 3.49 * ones(3,1);
        atom_angle = 90 * ones(3,1);
        atom_number = 3 * ones(2,1);
        [P, I, C, F] = deal(0, 1, 0, 0);
    case "Ba"
        atom_size = 5.02 * ones(3,1);
        atom_angle = 90 * ones(3,1);
        atom_number = 56 * ones(2,1);
        [P, I, C, F] = deal(0, 1, 0, 0);
    case "alpha-Fe"
        atom_size = 2.8664 * ones(3,1);
        atom_angle = 90 * ones(3,1);
        atom_number = 26 * ones(2,1);
        [P, I, C, F] = deal(0, 1, 0, 0);
    case "V"
        atom_size = 3.03 * ones(3,1);
        atom_angle = 90 * ones(3,1);
        atom_number = 23 * ones(2,1);
        [P, I, C, F] = deal(0, 1, 0, 0);
    case "Cr"
        atom_size = 2.8846 * ones(3,1);
        atom_angle = 90 * ones(3,1);
        atom_number = 24 * ones(2,1);
        [P, I, C, F] = deal(0, 1, 0, 0);
    case "Ta"
        atom_size = 3.30 * ones(3,1);
        atom_angle = 90 * ones(3,1);
        atom_number = 73 * ones(2,1);
        [P, I, C, F] = deal(0, 1, 0, 0);
    case "W"
        atom_size = 3.1650 * ones(3,1);
        atom_angle = 90 * ones(3,1);
        atom_number = 74 * ones(2,1);
        [P, I, C, F] = deal(0, 1, 0, 0);
    case "Cu"
        atom_size = 3.6147 * ones(3,1);
        atom_angle = 90 * ones(3,1);
        atom_number = 29 * ones(4,1);
        [P, I, C, F] = deal(0, 0, 0, 1);
    case "gamma-Fe"
        atom_size = 3.6468 * ones(3,1);%916 Centigrade
        atom_angle = 90 * ones(3,1);
        atom_number = 26 * ones(4,1);
        [P, I, C, F] = deal(0, 0, 0, 1);
    case "Ag"
        atom_size = 4.0857 * ones(3,1);
        atom_angle = 90 * ones(3,1);
        atom_number = 47 * ones(4,1);
        [P, I, C, F] = deal(0, 0, 0, 1);
    case "Al"
        atom_size = 4.0496 * ones(3,1);
        atom_angle = 90 * ones(3,1);
        atom_number = 13 * ones(4,1);
        [P, I, C, F] = deal(0, 0, 0, 1);
    case "Au"
        atom_size = 4.0788 * ones(3,1);
        atom_angle = 90 * ones(3,1);
        atom_number = 79 * ones(4,1);
        [P, I, C, F] = deal(0, 0, 0, 1);
    case "Ca"
        atom_size = 5.58 * ones(4,1);
        atom_angle = 90 * ones(3,1);
        atom_number = 20 * ones(4,1);
        [P, I, C, F] = deal(0, 0, 0, 1);
    case "Ni"
        atom_size = 3.5524 * ones(3,1);
        atom_angle = 90 * ones(3,1);
        atom_number = 28 * ones(4,1);
        [P, I, C, F] = deal(0, 0, 0, 1);
    otherwise
        error("材料数据库没有相关材料数据")
end
lattice_type = [P,I,C,F];
info = {atom_size, atom_angle, atom_number, lattice_type};
end

