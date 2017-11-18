function dateTimeStamp = getDateTimeStamp()
    dateTimeStamp = strrep(datestr(datetime('now')),'-','');
    dateTimeStamp = strrep(dateTimeStamp,' ','_');
    dateTimeStamp = strrep(dateTimeStamp,':','-');
end

