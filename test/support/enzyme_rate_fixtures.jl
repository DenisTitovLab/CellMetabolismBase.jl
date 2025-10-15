using CellMetabolismBase

function CellMetabolismBase.enzyme_rate(
    ::CellMetabolismBase.Enzyme{:Enz1,(:A_media,),(:A,)},
    metabs,
    params,
)
    return params.Enz1_Vmax * (metabs.A_media - metabs.A / params.Enz1_Keq) /
           (1 + metabs.A_media / params.Enz1_K_A_media + metabs.A / params.Enz1_K_A)
end

function CellMetabolismBase.enzyme_rate(
    ::CellMetabolismBase.Enzyme{:Enz2,(:A,),(:B, :B)},
    metabs,
    params,
)
    return params.Enz2_Vmax * (metabs.A - metabs.B^2 / params.Enz2_Keq) /
           (1 + metabs.A / params.Enz2_K_A + metabs.B^2 / params.Enz2_K_B)
end

function CellMetabolismBase.enzyme_rate(
    ::CellMetabolismBase.Enzyme{:Enz3,(:B,),(:C,)},
    metabs,
    params,
)
    return params.Enz3_Vmax * (metabs.B - metabs.C / params.Enz3_Keq) /
           (1 + metabs.B / params.Enz3_K_B + metabs.C / params.Enz3_K_C)
end

function CellMetabolismBase.enzyme_rate(
    ::CellMetabolismBase.Enzyme{:Enz4,(:C, :C),(:D,)},
    metabs,
    params,
)
    return params.Enz4_Vmax * (metabs.C^2 - metabs.D / params.Enz4_Keq) /
           (1 + metabs.C^2 / params.Enz4_K_C + metabs.D / params.Enz4_K_D)
end
