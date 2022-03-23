#merge several gating sets
data_dir <- "C:/Users/edmondsonef/Desktop/Humanized/Flow/"
wsp_file <- paste0(data_dir, '5-02Mar2022/15719 02Mar2022 Simone.wsp')
ws <- open_flowjo_xml(wsp_file)
ws


data_dir <- "C:/Users/edmondsonef/Desktop/Humanized/Flow/"
wsp_file1_05Jan2022 <- paste0(data_dir, '1-05Jan2022/15679 06Jan2022 Simone.wsp')
wsp_file2_02Feb2022 <- paste0(data_dir, '2-02Feb2022/15701 02Feb2022 Simone.wsp')
wsp_file3_15Feb2022 <- paste0(data_dir, '3-15Feb2022/15708 15Feb2022 Simone.wsp')
wsp_file4_28Feb2022 <- paste0(data_dir, '4-28Feb2022/15716 28Feb2022 Simone.wsp')
wsp_file5_02Mar2022 <- paste0(data_dir, '5-02Mar2022/15719 02Mar2022 Simone.wsp')


ws1 <- open_flowjo_xml(wsp_file1_05Jan2022)
ws1
ws2 <- open_flowjo_xml(wsp_file2_02Feb2022)
ws2
ws3 <- open_flowjo_xml(wsp_file3_15Feb2022)
ws3
ws4 <- open_flowjo_xml(wsp_file4_28Feb2022)
ws4

data_dir <- "C:/Users/edmondsonef/Desktop/Humanized/Flow/"
wsp_file5_02Mar2022 <- paste0(data_dir, '5-02Mar2022/15719 02Mar2022 Simone.wsp')
ws <- open_flowjo_xml(wsp_file5_02Mar2022)
ws
