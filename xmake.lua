add_rules("mode.debug", "mode.release")
set_config("mode", "debug")
set_config("arch", "x64")

add_requires("emhash","thread-pool")

function setup_target(target_name)
    set_kind("binary")
    set_languages("cxx20")
    set_rundir("$(projectdir)/")
    
    -- 平台特定的编译标志
    if is_plat("windows") then
        set_toolchains("msvc")
        add_cxflags("/utf-8")
    elseif is_plat("macosx", "linux") then
        set_toolchains("clang", "gcc")
        add_cxflags("-fPIC")
    end
    
    -- 调试模式下启用地址检查
    if is_mode("debug") then
        
    end
end

target("test_SA")
    setup_target("test_SA")
    add_files("optimizer/test_sa.cpp")
    add_packages("emhash")

target("test_BNB")
    setup_target("test_BNB")
    add_files("optimizer/test_bnb.cpp")
    add_packages("emhash")

target("test_TABU")
    setup_target("test_TABU")
    add_files("optimizer/test_tabu.cpp")
    add_packages("emhash")

target("test_MEM")
    setup_target("test_MEM")
    add_files("optimizer/test_mem.cpp")
    add_packages("emhash","thread-pool")

target("test_SA_MT")
    setup_target("test_SA_MT")
    add_files("optimizer/test_sa_mt.cpp")
    add_packages("emhash","thread-pool")